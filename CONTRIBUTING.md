# Contributing to ROED

Please open an issue before proposing a change to the statistical definition or
the optimization criterion. Code changes should preserve deterministic
single-threaded results and include tests that compare exact probabilities or
integer designs against independently verified values.

Run the package checks before submitting a pull request:

~~~r
R CMD build ROED
R CMD check --no-manual ROED_1.0.0.tar.gz
~~~

