function [] = play_run(mname)
rtwbuild(mname);
load_system(mname)
set_param(mname, 'SimulationCommand', 'connect')
set_param(mname, 'SimulationCommand', 'start')