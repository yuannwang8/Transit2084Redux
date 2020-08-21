# Transit2084Redux

This repo aims to visualise the 2084 transit of Earth and Moon from Mars from Martian orbit as well as various locations of interest on the surface.

Data courtesy of [NASA/JPL/HORIZONS](https://ssd.jpl.nasa.gov/?horizons).

This is the second attempt at this project: it improves upon the [attempt made in 2017](https://yuannwang8.github.io/Transit2084Mars/). Improvements so far:

* Now using the [`astropy`](https://www.astropy.org) package for common astrophysics tools to improve on calculation accuracies.
* Now using the [`astroquery`](https://astroquery.readthedocs.io/en/latest/) package to access the Horizons database. The list of interesting surface sites is curated in `MarsSurfaceSites.csv`.  Elevation details, where not provided, are estimated from [MOLA](https://attic.gsfc.nasa.gov/mola/images.html) imagery.

The sequence of activities can be found in `sequenceOfActivities.txt`. 

Example outputs are in the folder `ExampleOutputs`.
