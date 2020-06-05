# GeoGreensFunction

Commonly used Green's function in geoscience.

- [Line (2 nodes) antiplane dislocation]

- [Line (2 nodes) inplane dislocation]

- [Rectangular (4 nodes) dislocation (limited rotational freedom)]

- [Triangular (3 nodes) dislocation (singularity avoided)]

- [Rectangular (4 nodes) dislocation (fully rotational freedom)]

- [Rectangular (4 nodes) antiplane strain]

- [Rectangular (4 nodes) inplane strain]

- [Triangular (3 nodes) antiplane strain]

- [Triangular (3 nodes) inplane strain]

- [Hexahedron (8 nodes) volume strain]

- [Tetrahedron (4 nodes) volume strain]

All functions are translated from respective original sources (Matlab, Fortran, etc.), being validated,
pertaining all original parameter names and coordinate system. A transform to the ENU
(x -> east, y -> north, z -> upward) coordinate is provided if necessary. All methods here are type stable.
Please make a PR if you see anything is wrong or can be improved (^.^)!

*Some of them hasn't been implemented yet.*

## Reference

Okada, Y. (1992). Internal deformation due to shear and tensile faults in a half-space. Bulletin of the Seismological Society of America, 82(2), 1018–1040.

SEGALL, P. (2010). Earthquake and Volcano Deformation (STU-Student edition). Princeton University Press. Retrieved from http://www.jstor.org/stable/j.ctt7sg19

Nikkhoo, M., & Walter, T. R. (2015). Triangular dislocation: an analytical, artefact-free solution. Geophysical Journal International, 201(2), 1119–1141. https://doi.org/10.1093/gji/ggv035

Nikkhoo, M., Walter, T. R., Lundgren, P. R., & Prats-Iraola, P. (2017). Compound dislocation models (CDMs) for volcano deformation analyses. Geophysical Journal International, 208(2), 877–894. https://doi.org/10.1093/gji/ggw427

Barbot, S., Moore, J. D. P., & Lambert, V. (2017). Displacement and Stress Associated with Distributed Anelastic Deformation in a Half‐Space. Bulletin of the Seismological Society of America, 107(2), 821–855. https://doi.org/10.1785/0120160237

Barbot, S. (2018). Deformation of a Half‐Space from Anelastic Strain Confined in a Tetrahedral Volume. Bulletin of the Seismological Society of America, 108(5A), 2687–2712. https://doi.org/10.1785/0120180058
