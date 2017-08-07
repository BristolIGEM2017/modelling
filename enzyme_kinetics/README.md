# Enzyme Kinetics Folder
The enzyme kinetics model predicts the quantity of nitrite and nitrate reduced given their initial concentrations and the amound of Nap and Nrf enzymes present. The model uses the Michaelis-Menten enzyme-substrate chemical equation to create ODEs relating the concentrations of each species of the reaction over time. An example of the outputs of this model are given below (unrealistic concentrations and time scales):
![Michaelis-Menten Model](http://i.imgur.com/kAz6iM0.png)

The files in this folder are as follows:
* The **main** file `enzyme_kinetics.m` is supported by its ODE function file `enzyme_ODE.m`. This is the model used for the actual simulations, for example the image above.
* A **sensitivity study** on the rate constants can be performed with `run_sensitivity.m`, which is supported by a slightly modified version of the main file called `enzyme_kinetics_sensitivityfunction.m` (itself dependent on `enzyme_ode.m`).
* An **example** of Michaelis-Menten kinetics is given in `michaelis_menten_example.m`, which requires the ODE file `exampleODE.m`. It also compares the simulated ODEs with the quasi steady-state approximation (QSSA) to the Michaelis-Menten model, and therefore also requires the file `example_QSSA_ODE.m`.

Detail of the models will be given on the [Bristol iGEM wiki](http://2017.igem.org/Team:Bristol) at a later date.
