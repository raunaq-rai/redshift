"""
Interactive Redshift Fitting Tool
==================================
A simple tool to visually check emission line positions at different redshifts.
Use the slider to move where strong lines should appear on your spectrum.

Features:
- Main spectrum view with all emission line markers
- Zoom panels below for detailed view of each strong line region
- Fine-tuning slider with narrow range around z_prior
- Save button to record fitted redshifts to a file

Usage:
    from redshift_slider import RedshiftSlider
    
    # With your data
    slider = RedshiftSlider(wavelength, flux, z_prior=2.5, msaid='12345')
    slider.show()
    
    # With custom save file
    slider = RedshiftSlider(wavelength, flux, z_prior=2.5, msaid='12345',
                           save_file='my_redshifts.txt')
    slider.show()
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.gridspec import GridSpec
from pathlib import Path
from datetime import datetime

# Strong emission lines in rest-frame (Angstroms)
STRONG_LINES = {
    'Lyα': 1215.67,
    'NV': 1240.14,
    'CIV': 1549.06,
    'HeII': 1640.42,
    'CIII]': 1908.73,
    'MgII': 2798.75,
    '[OII]': 3727.09,
    'Hγ': 4340.47,
    'Hβ': 4861.33,
    '[OIII]4959': 4958.91,
    '[OIII]5007': 5006.84,
    'Hα': 6562.82,
    '[NII]6549': 6549.86,
    '[NII]6584': 6585.27,
    '[SII]6717': 6716.44,
    '[SII]6731': 6730.82,
}

# Key lines to show in zoom panels (most important for redshift fitting)
ZOOM_LINES = ['Lyα', 'CIV', 'CIII]', 'MgII', '[OII]', 'Hβ', '[OIII]5007', 'Hα']

# Colors for different line types
LINE_COLORS = {
    'Lyα': '#e74c3c',      # red
    'NV': '#e74c3c',
    'CIV': '#3498db',      # blue
    'HeII': '#3498db',
    'CIII]': '#3498db',
    'MgII': '#9b59b6',     # purple
    '[OII]': '#2ecc71',    # green
    'Hγ': '#f39c12',       # orange (Balmer)
    'Hβ': '#f39c12',
    '[OIII]4959': '#2ecc71',
    '[OIII]5007': '#2ecc71',
    'Hα': '#f39c12',
    '[NII]6549': '#1abc9c',  # teal
    '[NII]6584': '#1abc9c',
    '[SII]6717': '#95a5a6',  # gray
    '[SII]6731': '#95a5a6',
}


class RedshiftSlider:
    """
    Interactive redshift slider for visual line identification.
    
    Parameters
    ----------
    wavelength : array-like
        Observed wavelength array in Angstroms
    flux : array-like
        Flux array (any units)
    flux_err : array-like, optional
        Flux error array
    z_prior : float
        Initial redshift guess (default: 0.0)
    delta_z : float
        Slider range will be z_prior ± delta_z (default: 0.05)
    msaid : str, optional
        Object identifier (e.g., MSA ID). Used when saving redshifts.
    save_file : str or Path, optional
        Path to save fitted redshifts. Default: 'fitted_redshifts.txt'
    lines : dict, optional
        Custom line dictionary {name: rest_wavelength_A}
        If None, uses default STRONG_LINES
    zoom_lines : list, optional
        Which lines to show in zoom panels. Default: ZOOM_LINES
    zoom_width_A : float
        Width of zoom window in rest-frame Angstroms (default: 200)
    
    Attributes
    ----------
    z : float
        Current redshift value (updated when slider moves)
    saved : bool
        Whether the current redshift has been saved
    """
    
    def __init__(self, wavelength, flux, flux_err=None, z_prior=0.0,
                 delta_z=0.05, msaid=None, save_file='fitted_redshifts.txt',
                 lines=None, zoom_lines=None, zoom_width_A=200):
        
        self.wavelength = np.asarray(wavelength)
        self.flux = np.asarray(flux)
        self.flux_err = np.asarray(flux_err) if flux_err is not None else None
        self.z = z_prior
        self.z_prior = z_prior
        self.delta_z = delta_z
        self.zoom_width_A = zoom_width_A
        
        # Object ID and save file
        self.msaid = str(msaid) if msaid is not None else 'unknown'
        self.save_file = Path(save_file)
        self.saved = False
        
        # Set narrow z limits for fine tuning
        self.z_min = max(0.0, z_prior - delta_z)
        self.z_max = z_prior + delta_z
        
        # Use provided lines or default
        self.lines = lines if lines is not None else STRONG_LINES
        self.zoom_line_names = zoom_lines if zoom_lines is not None else ZOOM_LINES
        
        # Filter zoom lines to only those in wavelength range
        self._filter_zoom_lines()
        
        # Store line artists for updating
        self.line_artists = {}
        self.label_artists = {}
        self.zoom_line_artists = {}
        self.zoom_axes = {}
        
        self._setup_figure()
    
    def _filter_zoom_lines(self):
        """Keep only zoom lines that fall within the wavelength range."""
        wl_min, wl_max = self.wavelength.min(), self.wavelength.max()
        
        self.active_zoom_lines = []
        for name in self.zoom_line_names:
            if name in self.lines:
                rest_wl = self.lines[name]
                # Check if line is visible anywhere in the z range
                obs_wl_min = rest_wl * (1 + self.z_min)
                obs_wl_max = rest_wl * (1 + self.z_max)
                if obs_wl_max >= wl_min and obs_wl_min <= wl_max:
                    self.active_zoom_lines.append(name)
        
        # Limit to max 6 zoom panels
        self.active_zoom_lines = self.active_zoom_lines[:6]
    
    def _setup_figure(self):
        """Set up the matplotlib figure with slider and zoom panels."""
        n_zoom = len(self.active_zoom_lines)
        
        # Calculate figure layout
        if n_zoom == 0:
            # No zoom panels, simple layout
            self.fig, self.ax_main = plt.subplots(figsize=(14, 6))
            plt.subplots_adjust(bottom=0.15)
        else:
            # Create layout with zoom panels
            n_cols = min(n_zoom, 3)
            n_rows = (n_zoom + n_cols - 1) // n_cols
            
            self.fig = plt.figure(figsize=(14, 6 + 2.5 * n_rows))
            
            # GridSpec: main plot on top, zoom panels below, slider at bottom
            gs = GridSpec(2 + n_rows, n_cols, height_ratios=[3] + [1.5] * n_rows + [0.3],
                         hspace=0.4, wspace=0.3)
            
            # Main spectrum axis (spans all columns)
            self.ax_main = self.fig.add_subplot(gs[0, :])
            
            # Zoom panel axes
            for i, line_name in enumerate(self.active_zoom_lines):
                row = 1 + i // n_cols
                col = i % n_cols
                ax = self.fig.add_subplot(gs[row, col])
                self.zoom_axes[line_name] = ax
            
            plt.subplots_adjust(bottom=0.08, top=0.95, left=0.06, right=0.98)
        
        # Plot main spectrum
        self._plot_main_spectrum()
        
        # Plot zoom panels
        self._plot_zoom_panels()
        
        # Create slider
        ax_slider = plt.axes([0.12, 0.02, 0.45, 0.02])
        self.slider = Slider(
            ax=ax_slider,
            label=f'z (±{self.delta_z})',
            valmin=self.z_min,
            valmax=self.z_max,
            valinit=self.z,
            valstep=0.0001,
            color='#3498db'
        )
        self.slider.on_changed(self._update)
        
        # Reset button
        ax_reset = plt.axes([0.60, 0.015, 0.08, 0.025])
        self.reset_button = Button(ax_reset, 'Reset', color='lightgray')
        self.reset_button.on_clicked(self._reset)
        
        # Wider range button
        ax_wider = plt.axes([0.70, 0.015, 0.08, 0.025])
        self.wider_button = Button(ax_wider, 'Wider ±', color='lightgray')
        self.wider_button.on_clicked(self._widen_range)
        
        # Save button (green)
        ax_save = plt.axes([0.80, 0.015, 0.08, 0.025])
        self.save_button = Button(ax_save, 'Save z', color='#2ecc71', hovercolor='#27ae60')
        self.save_button.on_clicked(self._save_redshift)
        
        # Close button
        ax_close = plt.axes([0.90, 0.015, 0.08, 0.025])
        self.close_button = Button(ax_close, 'Done', color='#e74c3c', hovercolor='#c0392b')
        self.close_button.on_clicked(self._close)
    
    def _plot_main_spectrum(self):
        """Plot the main spectrum view."""
        ax = self.ax_main
        
        # Plot spectrum
        if self.flux_err is not None:
            ax.fill_between(
                self.wavelength,
                self.flux - self.flux_err,
                self.flux + self.flux_err,
                alpha=0.3, color='gray'
            )
        ax.plot(self.wavelength, self.flux, 'k-', lw=0.8)
        
        # Get y-limits
        valid = np.isfinite(self.flux)
        if valid.any():
            ymin, ymax = np.percentile(self.flux[valid], [1, 99])
            padding = (ymax - ymin) * 0.15
            self.ymin = ymin - padding
            self.ymax = ymax + padding * 1.5
        else:
            self.ymin, self.ymax = 0, 1
        
        ax.set_ylim(self.ymin, self.ymax)
        ax.set_xlim(self.wavelength.min(), self.wavelength.max())
        
        # Draw emission lines
        self._draw_main_lines()
        
        # Labels
        ax.set_xlabel('Observed Wavelength (Å)', fontsize=11)
        ax.set_ylabel('Flux', fontsize=11)
        
        # Title with MSAID
        title_text = f'MSAID: {self.msaid}  |  z = {self.z:.5f}'
        self.title = ax.set_title(title_text, fontsize=14, fontweight='bold')
        
        # Add legend
        self._add_legend()
    
    def _draw_main_lines(self):
        """Draw emission lines on main plot."""
        wl_min, wl_max = self.wavelength.min(), self.wavelength.max()
        
        for name, rest_wl in self.lines.items():
            obs_wl = rest_wl * (1 + self.z)
            
            if wl_min <= obs_wl <= wl_max:
                color = LINE_COLORS.get(name, '#7f8c8d')
                
                line = self.ax_main.axvline(
                    obs_wl, color=color, linestyle='--',
                    alpha=0.7, lw=1.2
                )
                self.line_artists[name] = line
                
                label = self.ax_main.text(
                    obs_wl, self.ymax * 0.95, name,
                    rotation=90, va='top', ha='right',
                    fontsize=7, color=color, fontweight='bold'
                )
                self.label_artists[name] = label
            else:
                self.line_artists[name] = None
                self.label_artists[name] = None
    
    def _plot_zoom_panels(self):
        """Plot the zoom panels around each key line."""
        for line_name, ax in self.zoom_axes.items():
            rest_wl = self.lines[line_name]
            obs_wl = rest_wl * (1 + self.z)
            
            # Zoom window in observed frame
            zoom_half_width = self.zoom_width_A * (1 + self.z) / 2
            wl_lo = obs_wl - zoom_half_width
            wl_hi = obs_wl + zoom_half_width
            
            # Get data in this range
            mask = (self.wavelength >= wl_lo) & (self.wavelength <= wl_hi)
            
            if mask.any():
                wl_zoom = self.wavelength[mask]
                fl_zoom = self.flux[mask]
                
                # Plot spectrum in zoom
                if self.flux_err is not None:
                    err_zoom = self.flux_err[mask]
                    ax.fill_between(wl_zoom, fl_zoom - err_zoom, fl_zoom + err_zoom,
                                   alpha=0.3, color='gray')
                ax.plot(wl_zoom, fl_zoom, 'k-', lw=1)
                
                # Auto y-limits based on data
                valid = np.isfinite(fl_zoom)
                if valid.any():
                    ymin, ymax = np.percentile(fl_zoom[valid], [2, 98])
                    padding = (ymax - ymin) * 0.2
                    ax.set_ylim(ymin - padding, ymax + padding)
            
            ax.set_xlim(wl_lo, wl_hi)
            
            # Draw the emission line
            color = LINE_COLORS.get(line_name, '#7f8c8d')
            zoom_line = ax.axvline(obs_wl, color=color, linestyle='-', lw=2, alpha=0.8)
            self.zoom_line_artists[line_name] = zoom_line
            
            # Title with line name
            ax.set_title(f'{line_name} ({rest_wl:.1f} Å)', fontsize=10, fontweight='bold',
                        color=color)
            ax.tick_params(labelsize=8)
    
    def _update(self, val):
        """Update line positions when slider changes."""
        self.z = self.slider.val
        self.saved = False  # Mark as unsaved when z changes
        wl_min, wl_max = self.wavelength.min(), self.wavelength.max()
        
        # Update main plot lines
        for name, rest_wl in self.lines.items():
            obs_wl = rest_wl * (1 + self.z)
            
            line = self.line_artists.get(name)
            label = self.label_artists.get(name)
            
            in_range = wl_min <= obs_wl <= wl_max
            
            if in_range:
                if line is None:
                    color = LINE_COLORS.get(name, '#7f8c8d')
                    line = self.ax_main.axvline(
                        obs_wl, color=color, linestyle='--',
                        alpha=0.7, lw=1.2
                    )
                    self.line_artists[name] = line
                    label = self.ax_main.text(
                        obs_wl, self.ymax * 0.95, name,
                        rotation=90, va='top', ha='right',
                        fontsize=7, color=color, fontweight='bold'
                    )
                    self.label_artists[name] = label
                else:
                    line.set_xdata([obs_wl, obs_wl])
                    line.set_visible(True)
                    label.set_x(obs_wl)
                    label.set_visible(True)
            else:
                if line is not None:
                    line.set_visible(False)
                if label is not None:
                    label.set_visible(False)
        
        # Update zoom panels
        for line_name, ax in self.zoom_axes.items():
            rest_wl = self.lines[line_name]
            obs_wl = rest_wl * (1 + self.z)
            
            # Update zoom window position
            zoom_half_width = self.zoom_width_A * (1 + self.z) / 2
            wl_lo = obs_wl - zoom_half_width
            wl_hi = obs_wl + zoom_half_width
            ax.set_xlim(wl_lo, wl_hi)
            
            # Update the vertical line position
            zoom_line = self.zoom_line_artists.get(line_name)
            if zoom_line:
                zoom_line.set_xdata([obs_wl, obs_wl])
            
            # Update y-limits based on new window
            mask = (self.wavelength >= wl_lo) & (self.wavelength <= wl_hi)
            if mask.any():
                fl_zoom = self.flux[mask]
                valid = np.isfinite(fl_zoom)
                if valid.any():
                    ymin, ymax = np.percentile(fl_zoom[valid], [2, 98])
                    padding = (ymax - ymin) * 0.2
                    ax.set_ylim(ymin - padding, ymax + padding)
        
        # Update title
        title_text = f'MSAID: {self.msaid}  |  z = {self.z:.5f}'
        self.title.set_text(title_text)
        self.fig.canvas.draw_idle()
    
    def _reset(self, event):
        """Reset slider to initial z_prior."""
        self.slider.set_val(self.z_prior)
    
    def _widen_range(self, event):
        """Widen the slider range."""
        self.delta_z *= 2
        self.z_min = max(0.0, self.z_prior - self.delta_z)
        self.z_max = self.z_prior + self.delta_z
        
        # Update slider
        self.slider.valmin = self.z_min
        self.slider.valmax = self.z_max
        self.slider.ax.set_xlim(self.z_min, self.z_max)
        self.slider.label.set_text(f'z (±{self.delta_z:.3f})')
        self.fig.canvas.draw_idle()
    
    def _save_redshift(self, event):
        """Save the current redshift to file."""
        # Create header if file doesn't exist
        if not self.save_file.exists():
            with open(self.save_file, 'w') as f:
                f.write("# Fitted Redshifts\n")
                f.write("# Generated by redshift_slider\n")
                f.write(f"# Created: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write("#\n")
                f.write("# msaid          z_prior       z_fitted      delta_z\n")
                f.write("#" + "-" * 60 + "\n")
        
        # Append the new redshift
        delta_z = self.z - self.z_prior
        with open(self.save_file, 'a') as f:
            f.write(f"{self.msaid:<16} {self.z_prior:<13.5f} {self.z:<13.5f} {delta_z:+.5f}\n")
        
        self.saved = True
        
        # Visual feedback - flash the save button
        self.save_button.color = '#27ae60'
        self.save_button.label.set_text('Saved!')
        self.fig.canvas.draw_idle()
        
        # Reset button after a moment
        self.fig.canvas.flush_events()
        import time
        time.sleep(0.5)
        self.save_button.color = '#2ecc71'
        self.save_button.label.set_text('Save z')
        self.fig.canvas.draw_idle()
        
        print(f"Saved: MSAID={self.msaid}, z={self.z:.5f} to {self.save_file}")
    
    def _close(self, event):
        """Close the figure."""
        if not self.saved and self.z != self.z_prior:
            print(f"Warning: Redshift changed but not saved! z={self.z:.5f}")
        plt.close(self.fig)
    
    def _add_legend(self):
        """Add a simplified legend for line types."""
        from matplotlib.lines import Line2D
        
        legend_elements = [
            Line2D([0], [0], color='#e74c3c', linestyle='--', label='UV (Lyα, NV)'),
            Line2D([0], [0], color='#3498db', linestyle='--', label='UV (CIV, CIII])'),
            Line2D([0], [0], color='#9b59b6', linestyle='--', label='MgII'),
            Line2D([0], [0], color='#2ecc71', linestyle='--', label='[OII], [OIII]'),
            Line2D([0], [0], color='#f39c12', linestyle='--', label='Balmer'),
        ]
        self.ax_main.legend(handles=legend_elements, loc='upper right', fontsize=8)
    
    def show(self):
        """Display the interactive plot."""
        plt.show()
        return self.z
    
    def get_redshift(self):
        """Return the current redshift value."""
        return self.z


def main():
    """Demo with synthetic data."""
    np.random.seed(42)
    
    # Observed wavelength range (simulating z~2 spectrum from NIRSpec)
    wavelength = np.linspace(10000, 25000, 5000)
    
    # Fake continuum + noise
    continuum = 1e-18 * (wavelength / 15000) ** (-1.5)
    noise = np.random.normal(0, 0.1 * continuum)
    
    # Add some fake emission lines at z=2.0
    z_true = 2.0
    for name, rest_wl in [('Hα', 6562.82), ('[OIII]5007', 5006.84), 
                          ('Hβ', 4861.33), ('[OII]', 3727.09)]:
        obs_wl = rest_wl * (1 + z_true)
        if wavelength.min() < obs_wl < wavelength.max():
            sigma = 50
            amplitude = 5e-18
            line = amplitude * np.exp(-0.5 * ((wavelength - obs_wl) / sigma) ** 2)
            continuum += line
    
    flux = continuum + noise
    flux_err = 0.1 * np.abs(continuum)
    
    print("=" * 60)
    print("Interactive Redshift Slider Demo")
    print("=" * 60)
    print(f"True redshift: z = {z_true:.4f}")
    print(f"Starting at z_prior = 1.98 (slightly off)")
    print("\nButtons:")
    print("  - Reset: Return to z_prior")
    print("  - Wider ±: Double the slider range")
    print("  - Save z: Save current redshift to file")
    print("  - Done: Close the window")
    print("=" * 60)
    
    slider = RedshiftSlider(
        wavelength, flux, flux_err=flux_err,
        z_prior=1.98,
        delta_z=0.05,
        msaid='demo_12345',
        save_file='fitted_redshifts.txt',
        zoom_width_A=150
    )
    final_z = slider.show()
    print(f"\nFinal redshift: z = {final_z:.5f}")
    print(f"Error from true: Δz = {final_z - z_true:.5f}")


if __name__ == '__main__':
    main()
