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
- Batch mode to process multiple spectra in sequence

Usage:
    from redshift_slider import RedshiftSlider, batch_fit
    
    # Single spectrum
    slider = RedshiftSlider(wavelength, flux, z_prior=2.5, msaid='12345')
    slider.show()
    
    # Batch mode - multiple spectra
    results = batch_fit(
        wavelengths=[wl1, wl2, wl3],
        fluxes=[fl1, fl2, fl3],
        z_priors=[2.5, 3.1, 1.8],
        msaids=['obj1', 'obj2', 'obj3']
    )
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
ZOOM_LINES = ['Lyα', 'CIV', 'CIII]', 'Hβ', '[OIII]5007', 'Hα']

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
        Width of zoom window in rest-frame Angstroms (default: 80, tighter zoom)
    figsize : tuple, optional
        Figure size (width, height). Default: (11, 5) for compact view
    batch_info : str, optional
        Text showing batch progress, e.g. "2/10"
    
    Attributes
    ----------
    z : float
        Current redshift value (updated when slider moves)
    saved : bool
        Whether the current redshift has been saved
    """
    
    def __init__(self, wavelength, flux, flux_err=None, z_prior=0.0,
                 delta_z=0.05, msaid=None, save_file='fitted_redshifts.txt',
                 lines=None, zoom_lines=None, zoom_width_A=80,
                 figsize=None, batch_info=None):
        
        self.wavelength = np.asarray(wavelength)
        self.flux = np.asarray(flux)
        self.flux_err = np.asarray(flux_err) if flux_err is not None else None
        self.z = z_prior
        self.z_prior = z_prior
        self.delta_z = delta_z
        self.zoom_width_A = zoom_width_A
        self.batch_info = batch_info
        
        # Figure size - smaller by default
        self.figsize = figsize if figsize is not None else (11, 5)
        
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
        
        # Calculate figure layout - compact sizing
        if n_zoom == 0:
            # No zoom panels, simple layout
            self.fig, self.ax_main = plt.subplots(figsize=self.figsize)
            plt.subplots_adjust(bottom=0.18)
        else:
            # Create layout with zoom panels
            n_cols = min(n_zoom, 3)
            n_rows = (n_zoom + n_cols - 1) // n_cols
            
            # Compact figure height
            fig_height = self.figsize[1] + 1.8 * n_rows
            self.fig = plt.figure(figsize=(self.figsize[0], fig_height))
            
            # GridSpec: main plot on top, zoom panels below, slider at bottom
            gs = GridSpec(2 + n_rows, n_cols, height_ratios=[2.5] + [1.2] * n_rows + [0.25],
                         hspace=0.45, wspace=0.3)
            
            # Main spectrum axis (spans all columns)
            self.ax_main = self.fig.add_subplot(gs[0, :])
            
            # Zoom panel axes
            for i, line_name in enumerate(self.active_zoom_lines):
                row = 1 + i // n_cols
                col = i % n_cols
                ax = self.fig.add_subplot(gs[row, col])
                self.zoom_axes[line_name] = ax
            
            plt.subplots_adjust(bottom=0.08, top=0.94, left=0.07, right=0.98)
        
        # Plot main spectrum
        self._plot_main_spectrum()
        
        # Plot zoom panels
        self._plot_zoom_panels()
        
        # Create slider
        ax_slider = plt.axes([0.12, 0.02, 0.40, 0.02])
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
        
        # Buttons - compact layout
        btn_y = 0.015
        btn_h = 0.025
        btn_w = 0.08
        
        # Reset button
        ax_reset = plt.axes([0.55, btn_y, btn_w, btn_h])
        self.reset_button = Button(ax_reset, 'Reset', color='lightgray')
        self.reset_button.on_clicked(self._reset)
        
        # Wider range button
        ax_wider = plt.axes([0.64, btn_y, btn_w, btn_h])
        self.wider_button = Button(ax_wider, 'Wider', color='lightgray')
        self.wider_button.on_clicked(self._widen_range)
        
        # Next/Done button (auto-saves) - green to indicate it saves
        ax_next = plt.axes([0.73, btn_y, btn_w + 0.02, btn_h])
        next_label = 'Next ✓' if self.batch_info else 'Done ✓'
        self.next_button = Button(ax_next, next_label, color='#2ecc71', hovercolor='#27ae60')
        self.next_button.on_clicked(self._next_and_save)
        
        # Skip button (orange) - skips without saving
        self.skipped = False
        if self.batch_info:
            ax_skip = plt.axes([0.84, btn_y, btn_w, btn_h])
            self.skip_button = Button(ax_skip, 'Skip', color='#e67e22', hovercolor='#d35400')
            self.skip_button.on_clicked(self._skip)
        else:
            self.skip_button = None
        
        # Connect keyboard events for arrow key navigation
        self.fig.canvas.mpl_connect('key_press_event', self._on_key_press)
    
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
        ax.plot(self.wavelength, self.flux, 'k-', lw=0.7)
        
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
        ax.set_xlabel('Observed Wavelength (Å)', fontsize=10)
        ax.set_ylabel('Flux', fontsize=10)
        ax.tick_params(labelsize=9)
        
        # Title with MSAID and batch info
        title_parts = []
        if self.batch_info:
            title_parts.append(f'[{self.batch_info}]')
        title_parts.append(f'MSAID: {self.msaid}')
        title_parts.append(f'z = {self.z:.5f}')
        
        title_text = '  |  '.join(title_parts)
        self.title = ax.set_title(title_text, fontsize=11, fontweight='bold')
        
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
                    alpha=0.7, lw=1.0
                )
                self.line_artists[name] = line
                
                label = self.ax_main.text(
                    obs_wl, self.ymax * 0.95, name,
                    rotation=90, va='top', ha='right',
                    fontsize=6, color=color, fontweight='bold'
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
            
            # Zoom window in observed frame - tighter zoom
            zoom_half_width = self.zoom_width_A * (1 + self.z) / 2
            wl_lo = obs_wl - zoom_half_width
            wl_hi = obs_wl + zoom_half_width
            
            # Get data in this range
            mask = (self.wavelength >= wl_lo) & (self.wavelength <= wl_hi)
            
            if mask.any():
                wl_zoom = self.wavelength[mask]
                fl_zoom = self.flux[mask]
                
                # Plot spectrum as steps (histogram style)
                if self.flux_err is not None:
                    err_zoom = self.flux_err[mask]
                    # Step-style error region
                    ax.fill_between(wl_zoom, fl_zoom - err_zoom, fl_zoom + err_zoom,
                                   alpha=0.3, color='gray', step='mid')
                ax.step(wl_zoom, fl_zoom, where='mid', color='k', lw=0.8)
                
                # Auto y-limits based on data - use full range to capture peaks
                valid = np.isfinite(fl_zoom)
                if valid.any():
                    ymin = np.nanmin(fl_zoom[valid])
                    ymax = np.nanmax(fl_zoom[valid])
                    padding = (ymax - ymin) * 0.15
                    ax.set_ylim(ymin - padding, ymax + padding)
            
            ax.set_xlim(wl_lo, wl_hi)
            
            # Draw the emission line - thinner
            color = LINE_COLORS.get(line_name, '#7f8c8d')
            zoom_line = ax.axvline(obs_wl, color=color, linestyle='-', lw=1, alpha=0.9)
            self.zoom_line_artists[line_name] = zoom_line
            
            # Title with line name
            ax.set_title(f'{line_name} ({rest_wl:.0f}Å)', fontsize=9, fontweight='bold',
                        color=color)
            ax.tick_params(labelsize=7)
    
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
                        alpha=0.7, lw=1.0
                    )
                    self.line_artists[name] = line
                    label = self.ax_main.text(
                        obs_wl, self.ymax * 0.95, name,
                        rotation=90, va='top', ha='right',
                        fontsize=6, color=color, fontweight='bold'
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
            
            # Update y-limits based on new window - use full range to capture peaks
            mask = (self.wavelength >= wl_lo) & (self.wavelength <= wl_hi)
            if mask.any():
                fl_zoom = self.flux[mask]
                valid = np.isfinite(fl_zoom)
                if valid.any():
                    ymin = np.nanmin(fl_zoom[valid])
                    ymax = np.nanmax(fl_zoom[valid])
                    padding = (ymax - ymin) * 0.15
                    ax.set_ylim(ymin - padding, ymax + padding)
        
        # Update title
        title_parts = []
        if self.batch_info:
            title_parts.append(f'[{self.batch_info}]')
        title_parts.append(f'MSAID: {self.msaid}')
        title_parts.append(f'z = {self.z:.5f}')
        self.title.set_text('  |  '.join(title_parts))
        self.fig.canvas.draw_idle()
    
    def _on_key_press(self, event):
        """Handle keyboard events for arrow key navigation."""
        if event.key == 'left':
            # Decrease z - regular step
            new_z = max(self.z_min, self.z - 0.0005)
            self.slider.set_val(new_z)
        elif event.key == 'right':
            # Increase z - regular step
            new_z = min(self.z_max, self.z + 0.0005)
            self.slider.set_val(new_z)
        elif event.key == 'shift+left':
            # Decrease z - fine step
            new_z = max(self.z_min, self.z - 0.0001)
            self.slider.set_val(new_z)
        elif event.key == 'shift+right':
            # Increase z - fine step
            new_z = min(self.z_max, self.z + 0.0001)
            self.slider.set_val(new_z)
        elif event.key == 'enter':
            # Save and move to next
            self._save_redshift()
            plt.close(self.fig)
    
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
    
    def _save_redshift(self):
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
        print(f"Saved: MSAID={self.msaid}, z={self.z:.5f} -> {self.save_file}")
    
    def _next_and_save(self, event):
        """Auto-save and close/move to next."""
        self._save_redshift()
        plt.close(self.fig)
    
    def _skip(self, event):
        """Skip this object without saving."""
        self.skipped = True
        print(f"Skipped: MSAID={self.msaid}")
        plt.close(self.fig)
    
    def _add_legend(self):
        """Add a simplified legend for line types."""
        from matplotlib.lines import Line2D
        
        legend_elements = [
            Line2D([0], [0], color='#e74c3c', linestyle='--', label='UV (Lyα)'),
            Line2D([0], [0], color='#3498db', linestyle='--', label='UV (CIV)'),
            Line2D([0], [0], color='#2ecc71', linestyle='--', label='[OIII]'),
            Line2D([0], [0], color='#f39c12', linestyle='--', label='Balmer'),
        ]
        self.ax_main.legend(handles=legend_elements, loc='upper right', fontsize=7)
    
    def show(self):
        """Display the interactive plot."""
        plt.show()
        return self.z
    
    def get_redshift(self):
        """Return the current redshift value."""
        return self.z


def batch_fit(wavelengths, fluxes, z_priors, msaids=None, flux_errs=None,
              save_file='fitted_redshifts.txt', **kwargs):
    """
    Process multiple spectra in sequence with the interactive slider.
    
    Parameters
    ----------
    wavelengths : list of arrays
        List of observed wavelength arrays
    fluxes : list of arrays
        List of flux arrays
    z_priors : list of float
        List of initial redshift guesses
    msaids : list of str, optional
        List of object identifiers. If None, uses indices.
    flux_errs : list of arrays, optional
        List of flux error arrays
    save_file : str or Path
        Path to save all fitted redshifts
    **kwargs
        Additional arguments passed to RedshiftSlider
        (e.g., delta_z, zoom_width_A, figsize)
    
    Returns
    -------
    results : dict
        Dictionary with keys:
        - 'msaid': list of object IDs
        - 'z_prior': list of prior redshifts
        - 'z_fitted': list of fitted redshifts
        - 'saved': list of booleans indicating if each was saved
        - 'skipped': list of booleans indicating if each was skipped
    
    Example
    -------
    >>> results = batch_fit(
    ...     wavelengths=[wl1, wl2, wl3],
    ...     fluxes=[fl1, fl2, fl3],
    ...     z_priors=[2.5, 3.1, 1.8],
    ...     msaids=['obj1', 'obj2', 'obj3'],
    ...     save_file='my_redshifts.txt'
    ... )
    """
    n_objects = len(wavelengths)
    
    # Validate inputs
    if len(fluxes) != n_objects:
        raise ValueError("wavelengths and fluxes must have same length")
    if len(z_priors) != n_objects:
        raise ValueError("wavelengths and z_priors must have same length")
    
    # Default MSAIDs
    if msaids is None:
        msaids = [f'obj_{i}' for i in range(n_objects)]
    elif len(msaids) != n_objects:
        raise ValueError("msaids must have same length as wavelengths")
    
    # Default flux errors
    if flux_errs is None:
        flux_errs = [None] * n_objects
    elif len(flux_errs) != n_objects:
        raise ValueError("flux_errs must have same length as wavelengths")
    
    # Results storage
    results = {
        'msaid': [],
        'z_prior': [],
        'z_fitted': [],
        'saved': [],
        'skipped': []
    }
    
    print("=" * 60)
    print(f"Batch Redshift Fitting: {n_objects} objects")
    print(f"Saving to: {save_file}")
    print("=" * 60)
    
    for i in range(n_objects):
        batch_info = f"{i+1}/{n_objects}"
        
        print(f"\n[{batch_info}] Processing MSAID: {msaids[i]} (z_prior={z_priors[i]:.4f})")
        
        slider = RedshiftSlider(
            wavelength=wavelengths[i],
            flux=fluxes[i],
            flux_err=flux_errs[i],
            z_prior=z_priors[i],
            msaid=msaids[i],
            save_file=save_file,
            batch_info=batch_info,
            **kwargs
        )
        
        z_fitted = slider.show()
        
        # Record results
        results['msaid'].append(msaids[i])
        results['z_prior'].append(z_priors[i])
        results['z_fitted'].append(z_fitted)
        results['saved'].append(slider.saved)
        results['skipped'].append(slider.skipped)
    
    # Summary
    n_saved = sum(results['saved'])
    n_skipped = sum(results['skipped'])
    print("\n" + "=" * 60)
    print(f"Batch Complete: {n_saved} saved, {n_skipped} skipped")
    print("=" * 60)
    
    return results


def main():
    """Demo with synthetic data - batch mode."""
    np.random.seed(42)
    
    # Generate 3 fake spectra at different redshifts
    z_trues = [2.0, 2.5, 1.8]
    wavelengths = []
    fluxes = []
    flux_errs = []
    z_priors = []
    msaids = []
    
    for i, z_true in enumerate(z_trues):
        # Observed wavelength range
        wl = np.linspace(8000, 26000, 5000)
        
        # Fake continuum + noise
        continuum = 1e-18 * (wl / 15000) ** (-1.5)
        noise = np.random.normal(0, 0.08 * continuum)
        
        # Add emission lines
        for name, rest_wl in [('Hα', 6562.82), ('[OIII]5007', 5006.84), 
                              ('Hβ', 4861.33), ('[OII]', 3727.09)]:
            obs_wl = rest_wl * (1 + z_true)
            if wl.min() < obs_wl < wl.max():
                sigma = 40
                amplitude = 4e-18
                line = amplitude * np.exp(-0.5 * ((wl - obs_wl) / sigma) ** 2)
                continuum += line
        
        wavelengths.append(wl)
        fluxes.append(continuum + noise)
        flux_errs.append(0.08 * np.abs(continuum))
        z_priors.append(z_true + np.random.uniform(-0.02, 0.02))  # Slightly off
        msaids.append(f'demo_{i+1}')
    
    print("=" * 60)
    print("Batch Redshift Slider Demo")
    print("=" * 60)
    print(f"True redshifts: {z_trues}")
    print("Starting with slightly offset z_priors")
    print("\nButtons:")
    print("  - Next ✓: Save and move to next object")
    print("  - Skip: Skip without saving")
    print("=" * 60)
    
    results = batch_fit(
        wavelengths=wavelengths,
        fluxes=fluxes,
        flux_errs=flux_errs,
        z_priors=z_priors,
        msaids=msaids,
        save_file='demo_redshifts.txt',
        delta_z=0.05,
        zoom_width_A=80
    )
    
    print("\nResults:")
    for i in range(len(results['msaid'])):
        status = "saved" if results['saved'][i] else ("skipped" if results['skipped'][i] else "not saved")
        print(f"  {results['msaid'][i]}: z={results['z_fitted'][i]:.5f} ({status})")


if __name__ == '__main__':
    main()
