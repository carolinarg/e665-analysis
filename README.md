# e665-analysis
#Access to BeAGLE:
If you have access to BNL servers you can directly use BeAGLE, see https://eic.github.io/software/beagle.html [1]
If you want to install BeAGLE locally, you have access from here: https://gitlab.in2p3.fr/BeAGLE/BeAGLE

Generation:
1. To set up your simulation for the E665 experiment, you need first set up your PYTHIA control card "S1ALL003", and you must set up your BeAGLE 
control card; it will be an input file, for example, "muXe_pyqm_qhat_g1.inp".
2. PyQM should be activated in the BeAGLE control card, as in the example "muXe_pyqm_qhat_g1.inp" or "muD_pyqm_qhat_g1.inp".
3. Following the instructions in [1], you generate the simulations, where you will obtain a text file that you need to convert to a root file.
4. With this root file, we are ready for the analysis!

Analysis:
1. To replicate the e665 multiplicity ratio analysis, you can use the "runEICTree_statis.C" for Xenon and "runEICTree_D_statis.C" for Deuterium.
2. Since for Deuterium, we only need the average number of particles, we create histograms for each rapidity region and charged particles, and we
take the average. We'll create a text file with these averages that you need to include in "runEICTree_statis.C", for example, 
"scale_g1_qhat_05_1E8.txt"[2].

3. We run the analysis for Xenon, after added the correspondly scales files[2] for each option, as 
      .L runEICTree_statis.C+
       runEICTree_statis.C()

4. (Old version) Obtaining the final ROOT file that we will use to overlap all the plots for each option, using "overlaping.C" and to get the final version
of the multiplicity ratios, we run "draw_final_plots.C".
    (New version) Plots are made using python.


