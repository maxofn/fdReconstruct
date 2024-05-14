# fdReconstruct

This repository provides the code and data considered in the work

Ofner, M., & Hörmann, S. (2024), "Covariate-informed reconstruction of partially observed functional data via factor models".

## Installation

The R-package `fdReconstruct` can be installed and loaded via the commands

```
devtools::install_github("maxofn/fdReconstruct/fdReconstruct")
library(fdReconstruct)
```

## Help

The main function for the reconstruction of functional data is `fdReconstruct`. Type

```
help(fdReconstruct)
```

for further details. Data of the real data application are contained in `temp_graz`.

## License

The package `fdReconstruct` is licensed under GNU General Public License v3.0.

## References

Kneip, A., & Liebl, D. (2020). On the optimal reconstruction of partially observed functional data. The Annals of Statistics, 48(3), 1692–1717.

Kraus, D. (2015). Components and completion of partially observed functional data. Journal of the Royal Statistical Society: Series B: Statistical Methodology, 777-801.

Land Steiermark (2023). Air quality data. https://app.luis.steiermark.at/luft2/suche.php, Licensed under CC BY 4.0; accessed on March 28, 2023.

Yao, F., Müller, H. G., & Wang, J. L. (2005). Functional data analysis for sparse longitudinal data. Journal of the American statistical association, 100(470), 577-590.
