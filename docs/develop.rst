
Develop
=======

Requires
--------

- `git <https://git-scm.com/>`__
- `Make <https://www.gnu.org/software/make/>`__. See `Makefile <https://github.com/contrailcirrus/pycontrails/blob/main/Makefile>`__ for a list of ``make`` commands.

Developing documentation requires:

- `pandoc <https://pandoc.org/installing.html>`__ for interpreting Jupyter notebooks
- `LaTeX <https://www.latex-project.org/get/>`__ for pdf outputs.
  If you are using a Mac, `MacTeX <https://www.tug.org/mactex/index.html>`__ is the best option.
  Note that LaTeX can be fairly large to install (~6GB).

Environment
-----------

Create a dedicated virtual environment for development:

.. code-block:: bash

    # create environment in <DIR>
    $ python3 -m venv <DIR>

    # activate environment (Unix-like)
    $ source <DIR>/bin/activate

If using `Anaconda <https://www.anaconda.com/>`__ / `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`__
Python, create a dedicated Anaconda environment:

.. code-block:: bash

    # create conda environment
    $ conda create -n contrails python=3.10

    # activate environment
    $ conda activate contrails


Install
-------

After activating the virtual environment, clone the `pycontrails repository <https://github.com/contrailcirrus/pycontrails>`__:

.. code-block:: bash

    $ cd <install-path>
    $ git clone git@github.com:contrailcirrus/pycontrails.git
    $ cd pycontrails

Install the development verison of ``pycontrails`` using ``make``:

.. code-block:: bash

    $ make dev-install

Install dependencies manually using ``pip`` in editable mode:

.. code-block:: bash

    # core development installation
    $ pip install -e ".[docs,dev]"

    # install optional dependencies as above
    $ pip install -e ".[ecmwf,gfs]"

    # make sure to add the pre-commit hooks if installing manually
    $ pre-commit install


Test
----

Run all code quality checks and unit tests:

.. code-block:: bash

    $ make test

Lint the repository with ``flake8``:

.. code-block:: bash

    $ make flake8

Autoformat the repository with ``black``:

.. code-block:: bash

    $ make black

Run type checking with ``mypy``:

.. code-block:: bash

    $ make mypy

Run unit tests with ``pytest``:

.. code-block:: bash

    $ make pytest

Documentation
-------------

Documentation is written in `reStructuredText <http://docutils.sourceforge.net/rst.html>`__
and built with `Sphinx <https://www.sphinx-doc.org/en/master/>`__.

Sphinx configuration is written in `docs/conf.py <https://github.com/contrailcirrus/pycontrails/blob/main/docs/conf.py>`__.
See `Sphinx configuration docs <https://www.sphinx-doc.org/en/master/usage/configuration.html>`__ for the full list of configuration options.

Build HTML documentation:

.. code-block:: bash

    # docs build to directory docs/_build/html
    $ make docs-build

    # automatically build docs on changes
    # docs will be served at http://127.0.0.1:8000
    $ make docs-serve

    # cleanup all built documentation
    $ make docs-clean

Build manually with ``sphinx-build``:

.. code-block:: bash

    $ sphinx-build -b html docs docs/_build/html      # HTML output

Sphinx caches builds between changes.
To force the whole site to rebuild, use the options ``-aE``:

.. code-block:: bash

    $ sphinx-build -aE -b html docs docs/_build/html  # rebuild all output

See `sphinx-build <https://www.sphinx-doc.org/en/master/man/sphinx-build.html#cmdoption-sphinx-build-b>`__
for a list of all the possible output builders.

PDF Output
~~~~~~~~~~

    Building PDF output requires a `LaTeX distribution <https://www.latex-project.org/get/>`__.

Build pdf documentation:

.. code-block:: bash

    $ make docs-pdf

A single pdf output (i.e. ``pycontrails.pdf``) will be built within ``docs/_build/latex``.

To build manually, run:


.. code-block:: bash

    $ sphinx-build -b latex docs docs/_build/latex
    $ cd docs/_build/latex
    $ make

References
~~~~~~~~~~

Bibliography references managed in a `Zotero library <https://www.zotero.org/groups/4730892/pycontrails/library>`__.

To automatically sync this library with the
`docs/_static/pycontrails.bib <https://github.com/contrailcirrus/pycontrails/blob/main/docs/_static/pycontrails.bib>`__ Bibtex file:

- Install `Zotero <https://www.zotero.org/>`__ and add the `pycontrails collection <https://www.zotero.org/groups/4730892/pycontrails/library>`__.
- Install the `Zotero Better Bibtex extension <https://retorque.re/zotero-better-bibtex/installation/>`__. Leave defaults during setup.
- Right click on the **pycontrails** library and select *Export Library*
- Export as *Better Bibtex*. You can optionally check *Keep Updated* if you want
  this file to update every time you make a change to the library.
- Select the file ``_static/pycontrails.bib`` and press *Save* to overwrite the file.
- Commit the updated ``_static/pycontrails.bib``

Test
~~~~

    All doc tests first ensure ERA5 data is cached locally:

    .. code-block:: bash

        $ make ensure-era5-cached

Run docstring example tests with `doctest <https://docs.python.org/3/library/doctest.html>`__:

.. code-block:: bash

    $ make doctest

Test notebook examples with `nbval pytest plugin <https://github.com/computationalmodelling/nbval>`__:

.. code:: bash

   $ make nbtest
