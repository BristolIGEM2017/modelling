# Bristol iGEM 2017 Modelling
This repository contains the mathematical models and simulations for the single-cell gene expression model and enzyme kinetics model, as well as the multi-cell population model using Bristol iGEM 2008's BSim tool. Each aspect is grouped into its respective folder, which contains the relevant readme.

## Gene Expression
This model takes IPTG (inducer) concentration as an input and outputs the concentrations of Nap and Nrf enzymes produced. Modelling gene expression using ordinary differential equations (ODEs) allows for time-accurate representations of enzyme production given impulses of IPTG.

## Enzyme Kinetics
The enzyme kinetics model predicts the quantity of nitrite and nitrate reduced given their initial concentrations and the amound of Nap and Nrf enzymes present. The model uses the Michaelis-Menten enzyme-substrate chemical equation to create ODEs relating the concentrations of each species of the reaction over time.
![Michaelis-Menten Model](http://i.imgur.com/kAz6iM0.png)

## Population Modelling
BSim, an agent-based tool developed by Bristol's 2008 iGEM team, allows for propagation in time of bacteria, which can be defined to a high level of realism thanks to in-built ODE solvers and chemical field propagators amongst other features. Most importantly, it is simple to expand to fit the needs of any particular model, such as the implementation of the above ODEs.
