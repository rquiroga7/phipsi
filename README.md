# PhiPsi Viewer — Ramachandran Principle Clone

Static site cloning the 3D graphics and controls from **https://bioinformatics.org/molvis/phipsi/** for hosting on **GitHub Pages**.

- **3D viewer**: 3Dmol.js rendering the same tripeptide `phipsi-16atoms.pdb` (Lys15-Ala16-Arg17, 16 atoms, central Ala16). Same colors as original: C `#c8c8c8`, Cα `#707070`, N `#6580FF`, O `#FF6060`, bonds white, Phi bond green `#80ff80`, peptide bonds magenta `#ff80ff` when enabled, background `#d0d0d0`.
- **Same controls** as original (custom checkbox styling from `input.css`): **Alanine** (black dots 60%), **Peptide Bonds** (magenta), **Planes** (peptide planes), **van der Waals** (88% translucent, **White** toggle), **Show Clashes** (orange contacts, **Trail Clashes**), **Reset**. Same button styling `.rotationbutton` `#d8ffd8`/`#4CAF50`.
- **New Ramachandran heatplot**: Plotly heatmap (alpha/beta/left-handed Gaussians) with **red dot** and **red X/Y crosshairs** tracking current Phi/Psi. **Click anywhere** on the plot to animate the structure gradually to that Phi/Psi (same gradual animation as the +/- buttons).
- **Four buttons** (instead of original Phi/Psi selector + +/-20°): `+15 Phi`, `-15 Phi`, `+15 Psi`, `-15 Psi` — all animate in 2° steps (10 steps) with live Phi/Psi display and plot marker update. Non-3D tutorial text and static JPG plot removed; only 3D viewer + controls + interactive plot remain.

## Local use

Open `index.html` via a local http server (required for `fetch('phipsi-16atoms.pdb')`; `file://` also works via embedded fallback):

```sh
npx http-server . -p 8080
# then open http://localhost:8080
```

Or just open `index.html` directly — embedded PDB fallback ensures it still loads.

## GitHub Pages

This repo is ready for Pages:

1. Push `main` to `github.com/rquiroga7/phipsi`
2. In **Settings → Pages**: Source = `Deploy from a branch`, Branch = `main` / `/(root)`, Save.
3. Add `.nojekyll` (already included) so `_`-files and PDB are served.
4. Visit `https://rquiroga7.github.io/phipsi/`.

Original tutorial by Eric Martz (CC BY-NC-SA 4.0), JSmol by Bob Hanson. This clone uses 3Dmol.js + Plotly for a lightweight GitHub Pages deployment. See `phipsi-16atoms.pdb` (original) and `Ramachandran_plot_general_100K.jpg` (original static plot, now replaced by interactive heatplot).
