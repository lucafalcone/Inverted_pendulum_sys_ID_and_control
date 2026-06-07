function [] = play_run(mname)
mname = 'inverted_pendulum_template';
rtwbuild(mname);
load_system(mname)
set_param(mname, 'SimulationCommand', 'connect')
set_param(mname, 'SimulationCommand', 'start')