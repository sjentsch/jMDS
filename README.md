# jMDS: A jamovi Module for Multidimensional Scaling
jMDS is a jamovi-module for multidimensional scaling. It is available from the
[jamovi library](https://docs.jamovi.org/_pages/Install_modules.html).

## Install in jamovi

If you haven't already installed [jamovi](https://www.jamovi.org/download.html),
download and install it. Open jamovi, and select the jamovi modules library and
[install](https://jamovi.readthedocs.io/en/latest/howto/howto_Install_modules.html)
jMDS from there.

## Install the development version (from source)

Assuming that you have [jamovi](https://www.jamovi.org/download.html) and 
[R](https://cran.r-project.org/) already installed on your machine, open
`R` and install the R-package `jmvtools` using the following command (this
is only required once before the first use):

```
install.packages('jmvtools', repos='https://repo.jamovi.org')
```


Clone this repository (or [download the ZIP file](../../archive/refs/heads/main.zip)
and extract it). In `R`, go to the directory where you cloned or extracted the
repository to, and compile the module with:

```
jmvtools::install()
```

The module then appears among your jamovi modules.


## Giving back

If you find this module helpful, please consider donating to the jamovi project (via
the Patreon-link on the left side).
If you can't give money, but would like to support us in another way, you may contribute
to translating [jamovi](https://hosted.weblate.org/engage/jamovi/), the
[jamovi documentation](https://hosted.weblate.org/engage/jamovidocs/) or the textbook
[”learning statistics with jamovi“](https://hosted.weblate.org/engage/jamovi/) into your
language.

Thank you for your support!

