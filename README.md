# redshift_slider

Interactive redshift fitting tool with visual emission line markers.

## Installation

```bash
# From the package directory
pip install -e .
```

## Quick Start

```python
from redshift_slider import RedshiftSlider

# Load your spectrum (observed wavelength in Angstroms)
# wavelength, flux = ... your data ...

# Create the slider
slider = RedshiftSlider(
    wavelength, 
    flux,
    flux_err=flux_err,   # optional
    z_prior=2.5,         # your initial redshift guess
    msaid='12345',       # object identifier
    save_file='fitted_redshifts.txt'  # where to save results
)

# Show the interactive plot
z_fitted = slider.show()
```

## Features

- **Interactive slider** to adjust redshift and see where emission lines would fall
- **Zoom panels** around key emission lines for detailed inspection
- **Save button** to record fitted redshifts to a file with MSAID
- **Color-coded lines** by type (UV, Balmer, [OIII], etc.)
- Fine-tuning with **narrow z range** (±0.05 default) and 0.0001 step size

## Buttons

| Button | Function |
|--------|----------|
| Reset | Return to z_prior |
| Wider ± | Double the slider range |
| Save z | Save current redshift to file |
| Done | Close the window |

## Output File Format

The save file (`fitted_redshifts.txt`) contains:
```
# msaid          z_prior       z_fitted      delta_z
demo_12345       1.98000       2.00123       +0.02123
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `wavelength` | required | Observed wavelength array (Å) |
| `flux` | required | Flux array |
| `flux_err` | None | Flux error array |
| `z_prior` | 0.0 | Initial redshift guess |
| `delta_z` | 0.05 | Slider range: z_prior ± delta_z |
| `msaid` | 'unknown' | Object identifier |
| `save_file` | 'fitted_redshifts.txt' | Output file path |
| `zoom_width_A` | 200 | Zoom panel width (rest-frame Å) |

## Custom Lines

```python
# Use only specific lines
my_lines = {
    'Hα': 6562.82,
    '[OIII]5007': 5006.84,
    'Hβ': 4861.33,
}
slider = RedshiftSlider(wavelength, flux, z_prior=2.0, lines=my_lines)
```

## Demo

Run the built-in demo:
```bash
redshift-slider-demo
```

Or in Python:
```python
from redshift_slider.core import main
main()
```
