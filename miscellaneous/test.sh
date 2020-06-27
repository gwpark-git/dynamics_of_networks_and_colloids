# ./Brownian_simulation_one test.inp
python post_processing/plot_association_2d_onestep.py data_longer/test data_longer/figures_test 40 1000 12 10.0
ffmpeg -r 10 -i data_longer/figures_test/t%08d.png -vcodec h264 -vf scale=2400:-1 -pix_fmt yuv420p data_longer/test.mp4
