TGCA subtype classification
==============================
> Detecting cancer subtypes with machine learning.

This repository contains the data, code, and manuscript accompanying the
preprint:

> WF Flynn, S Namburi, CA Paisie, HV Reddi, S Li, KRK Murthy, J Georgy.
> "Pan-cancer machine learning predictors of tissue of origin and molecular
> subtype." Submitted, 2018.

currently available at
[bioRxiv](https://www.biorxiv.org/content/early/2018/05/30/333914).

## License
The code present in this repository is free to use for academic and
non-commercial use, and is subject to the following [License](LICENSE) (also
available in [`docx` format](LICENSE.docs).

## Project Organization
This project is organized using a subset of the [Cookiecutter Data
Science](https://drivendata.github.io/cookiecutter-data-science/#directory-structure)
project structure.

All data and results, and most visualizations can be generated from scratch
using the `make` command.  A full build of the project can be done with

```bash
make requirements
make data
make models
make viz
```

## Requirements
In order to produce the models and visualizations, this project requires
`conda`, through which `R` and `Python 3.6` will be installed along with their
needed modules/packages.

Running `make requirements` will:

* Create and activate a conda environment named `tcga_subtype_classification`.
* Install `R` and `Python 3.6` along with the packages listed in
  `requirements.txt` and `requirements_conda.txt`.
* Test these installations.

If you do not have a `conda` installation, you can install a minimal
installation through [`miniconda`](https://conda.io/miniconda.html).

## Figures and web application
Figures present in the manuscript preprint can be generated automatically using
`make viz` or interactively using the notebooks (symlinked) in the `/notebooks/`
root directory.

We've also include a simple interactive web vizualization that is currently
hosted at [Pan Cancer Classification Portal](https://pccportal.jax.org).  You
can host your own version locally using the code in the `/app/` root directory:

```bash
cd app/
python3 run_flask.py [--host IP] [--port PORT]
```

<p><small>Project based on the <a target="_blank"
href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data
science project template</a>. #cookiecutterdatascience</small></p>
