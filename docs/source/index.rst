.. torch-adata documentation master file.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: ../imgs/torch-adata.logo.large.svg
  :width: 1200
  :alt: torch-adata logo
======================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Documentation
********************
``torch-adata`` is a framework for bridging data held in `AnnData` - the most popular single-cell python-based data structure and companion to scanpy - to the PyTorch ``Dataset`` class. ``torch-adata`` is meant to be structured but flexible within the rules of both of these data structures while also being easy to use.

* Deep learning enthusiasts: get started faster with single-cell data.* Biologists: get started faster with ML
* Store and access single-cell data efficiently and share between models, reproducibly (echoing the goals of the ``LightningDataModule``.


Install torch-adata
-------------------

.. raw:: html

   <div class="row" style='font-size: 16px'>
      <div class='col-md-6'>

To install the PYPI distribution:

.. code-block:: bash

    pip install torch-adata

.. raw:: html

      </div>
      <div class='col-md-6'>

to install the development version:

.. code-block:: bash

    git clone https://github.com/mvinyard/torch-adata.git; cd torch-adata
    pip install -e .

.. raw:: html

      </div>
   </div>

.. toctree::
   :maxdepth: 3
   :titlesonly:
   :template: classtemplate.rst
   :hidden:

   installation
   API/index
   API/ancilliary_functions

.. note::

   This project is under active development.
   
A bit about me
--------------
I'm currently working on my PhD at Harvard University where I am developing methods for the analysis of single-cell genomics data. I also work on methods for CRISPR-Cas9-based genome-editing tools (like base-editing) to study cancer progression. My PhD work is funded by an NIH F31 award through the NCI and focuses on applying innovations in deep learning to analogous problems in the analysis of single-cell data.


If you find this package useful for your own work, please consider reaching out and connecting through `GitHub <https://github.com/mvinyard>`_, `LinkedIn <https://www.linkedin.com/in/michaelvinyard/>`_, or `Twitter <https://twitter.com/vinyard_m>`_. 

