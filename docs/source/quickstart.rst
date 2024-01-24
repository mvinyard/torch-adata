===========
Quick Start
===========

Install the library
-------------------

.. code-block:: shell

   pip install neural-diffeqs


Import the library
------------------

.. code-block:: python

   import neural_diffeqs


Data
----

Let's say you have some temporally-resolved data, spanning three time points. There are 200 12-dimension samples at each time point.

.. code-block:: python

   # the initial state
   X0 = torch.randn([200, 12])
   
   # later states
   X1 = torch.randn([200, 12])
   X2 = torch.randn([200, 12])
   
   time = torch.Tensor([1, 4, 7])
   

Neural network function
-----------------------

Call together an SDE (or ODE):

.. code-block:: python

   SDE = neural_diffeqs.NeuralSDE(
       state_size = 12, mu_hidden = [32, 32], sigma_hidden = [32, 32],
   )
   
   ODE = neural_diffeqs.NeuralSDE(state_size = 12, mu_hidden = [32, 32])


Make a prediction
-----------------

.. code-block:: python
   
   import torchsde
   
   # dt is an important parameter to tune; we'll start with 0.1
   
   X_hat_sde = torchsde.sdeint(SDE, X0, ts = time, dt = 0.1)
   
   X_hat_ode = torchsde.sdeint(ODE, X0, ts = time, dt = 0.1)
   
   # compare to X1, X2 (observed states)
   # integrate into a training regimen to fit the mu, sigma networks
