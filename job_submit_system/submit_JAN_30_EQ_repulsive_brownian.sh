N_runs=10
init_sample_number=0
bash automation_submit.sh equilibrium_detail/repulsive_brownian/NP0400_LD10P3_C025.inp $N_runs $init_sample_number ../LD10P3VER/EQ_NP0400_LD10P3/NP0400_LD10P3.traj
bash automation_submit.sh equilibrium_detail/repulsive_brownian/NP0600_LD10P3_C025.inp $N_runs $init_sample_number ../LD10P3VER/EQ_NP0600_LD10P3/NP0600_LD10P3.traj
bash automation_submit.sh equilibrium_detail/repulsive_brownian/NP0512_LD08P3_C025.inp $N_runs $init_sample_number ../LD08P3VER/EQ_NP0512_LD08P3/NP0512_LD08P3.traj

bash automation_submit.sh equilibrium_detail/repulsive_brownian/NP0400_LD10P3_C100.inp $N_runs $init_sample_number ../LD10P3VER/EQ_NP0400_LD10P3/NP0400_LD10P3.traj
bash automation_submit.sh equilibrium_detail/repulsive_brownian/NP0600_LD10P3_C100.inp $N_runs $init_sample_number ../LD10P3VER/EQ_NP0600_LD10P3/NP0600_LD10P3.traj
bash automation_submit.sh equilibrium_detail/repulsive_brownian/NP0512_LD08P3_C100.inp $N_runs $init_sample_number ../LD08P3VER/EQ_NP0512_LD08P3/NP0512_LD08P3.traj
