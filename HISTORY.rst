=======
History
=======

0.1.0 (2021-03-18)
------------------


0.1.3 (2021-09-31)
------------------
new implementation of case-creation.

-> try out an example with case-creation. create geometry and then case
-> put an ascii-based simulation into the template-dir
-> call the simulation by name of the directory in your configuration-file
-> use <var VAR_VAR var> and <opt OPT_OPT opt> as variable/option signature
-> define variables and options in your configuration-file


implementation of new testfunction-module.

->it is running all example-cases as defined and it takes time

partially fixed example-cases. need new version-incrementation for a finalized set of examples

0.1.4 (2021-10-04)
------------------
implementation of parastudy-creation. it is still necessary to implement options for the case-creation

0.1.5 (2021-10-14)
------------------
new feature: job-management. you can choose from a set of job-scripts. You can easily set parameters in your job-script
This allowes a simulation management for the parameter-study-method.

0.1.10 (2021-11-25)
------------------
loads of major updates.
    -new extract_hk_vk: way more reliable, easier to maintain and faster
    -implementation of a new profile-generator
    -generic postprocess caller-method
    -advanced job-management
