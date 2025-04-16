Installation
============

You can install the project dependencies using either **pip** (via `requirements.txt`) or **conda** (via `environment.yml`), depending on your preferred environment manager.

.. contents::
   :local:
   :depth: 2


Installation using pip
-----------------------

To install the required dependencies with pip:

1. (Optional) Create and activate a virtual environment:

   .. code-block:: bash

      python -m venv venv
      source venv/bin/activate 


2. Install the dependencies:

   .. code-block:: bash

      pip install -r docs/requirements.txt



Installation using conda
-------------------------

To set up the environment using conda:


1. Create the environment from the `environment.yml` file:

   .. code-block:: bash

      conda env create -f docs/environment.yml


2. Activate the environment:

   .. code-block:: bash

      conda activate wavepropenv


Verifying installation
----------------------

After installation, verify that everything is set up correctly by running the following:

1. Open Python:

   .. code-block:: bash

      python

2. Run the following command to ensure the package and functions are available:

   .. code-block:: python

      from scripts import *

   If no errors occur, the installation was successful, and you can begin using the functions from the **scripts** package to run simulations and plot model results. For examples, refer to :doc:`examples`.
