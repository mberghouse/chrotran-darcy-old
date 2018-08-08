#!/bin/bash

pflotran=../../../../../src/pflotran/pflotran

cd reaction_only/analytical
pwd
$pflotran -input_prefix decay_only_analytical >& decay_only_analytical.stdout
cd ../conventional
pwd
$pflotran -input_prefix decay_only_conventional >& decay_only_conventional.stdout
cd ../conventional_sm_ts
pwd
$pflotran -input_prefix decay_only_conventional >& decay_only_conventional.stdout
cd ../implicit
pwd
$pflotran -input_prefix decay_only_implicit >& decay_only_implicit.stdout
#------------------------------------------------------------------------------#
cd ../../reactive_transport/1d_structured/
pwd
$pflotran -input_prefix rt >& rt.stdout
cd ../1d_structured_uniform/
pwd
$pflotran -input_prefix rt >& rt.stdout
cd ../1d_unstructured_implicit/
pwd
$pflotran -input_prefix rt >& rt.stdout
cd ../1d_unstructured_irregular/
pwd
$pflotran -input_prefix rt >& rt.stdout
cd ../1d_unstructured_test/
pwd
$pflotran -input_prefix rt >& rt.stdout
cd ../1d_unstructured_uniform/
pwd
$pflotran -input_prefix rt >& rt.stdout
#------------------------------------------------------------------------------#
cd ../../transport_only/1d_structured/
pwd
$pflotran -input_prefix tracer >& tracer.stdout
cd ../1d_structured_uniform/
pwd
$pflotran -input_prefix tracer >& tracer.stdout
cd ../1d_unstructured_implicit/
pwd
$pflotran -input_prefix tracer >& tracer.stdout
cd ../1d_unstructured_irregular/
pwd
$pflotran -input_prefix tracer >& tracer.stdout
cd ../1d_unstructured_test/
pwd
$pflotran -input_prefix tracer >& tracer.stdout
cd ../1d_unstructured_uniform/
pwd
$pflotran -input_prefix tracer >& tracer.stdout
cd ../..
find . | grep stdout
