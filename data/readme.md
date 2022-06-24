# Symbolic Regression example datasets

These are the example datasets for use with "v1.0.0_YourSymbolicRegression-Master" notebook in the geppy examples directory.

## Example Symbolic Regression datasets

### UCI Power Plant
https://archive.ics.uci.edu/ml/datasets/combined+cycle+power+plant

### UCI Concrete Compressive Strength Data Set
https://archive.ics.uci.edu/ml/datasets/concrete+compressive+strength

This data problem proved to be hard for symbolic regression, as the equations that lower mse generate nan, inf, -inf, and zoo in many holdout cases.

It is provided to help illustrate these difficulties.
You will notice that there are options to work with z-score normalised data in the notebook. This was one of the ideas I tried to see if it could help.

### Nasa Airfoil Self-Noise Data Set 
https://archive.ics.uci.edu/ml/datasets/Airfoil+Self-Noise

Good notes about the data and performance of regression models is here:

https://claudinei-daitx.github.io/airfoil-self-noise-prediction/#about

### Notes

These three datasets can be all used with the symbolic regression notebook called "v1.0.0_YourSymbolicRegression-Master" if you configure it to load both the data and the dictionary files included.

If you have your own dataset, prepare it similarly to these examples, and provide a data dictionary, and the notebook should be easy to use for your new dataset.
- Use headers in your delimited data
- Describe the data in a dictionary file, as seen here.
- the dictionary file drives the automated configuration of the notebook! 

