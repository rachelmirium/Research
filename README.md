Preliminary code for a speech synthesizer from electromagnetic articulography or MRI data, based on work by Maeda in 1982 and 1996, and Toutios and Maeda in 2012.

VCV_nasal_inputs takes input files from a specified folder corresponding to a vowel-consonant-vowel combination, and a sampling frequency (use 10000 for best results)

synth_with_nt uses output from VCV_nasal_inputs to create scaled speech signal

Example to produce the sound /asi/:
[A, X, Ag0, Ap, F0, Anc, AN, XN] = VCV_nasal_inputs(asi-area, 10000);
signal = synth_with_nt(A, X, Ag0, Ap, F0, Anc, AN, XN, Fs);
