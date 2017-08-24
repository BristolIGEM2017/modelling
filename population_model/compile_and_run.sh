#!/bin/bash

javac -cp .:../lib/core.jar:../lib/vecmath.jar:../lib/objimport.jar:../lib/bsim.jar BSimEColi.java
java  -cp ..:../lib/core.jar:../lib/vecmath.jar:../lib/objimport.jar:../lib/bsim.jar population_model.BSimEColi
rm *.class
