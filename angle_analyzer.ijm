input = "\\\\WURNET.NL\\Homes\\gobes001\\My Documents\\To send\\test\\";
setBatchMode(true);
list = getFileList(input);
for (i = 0; i < list.length; i++){
		open(input + list[i]);
		if (nImages>=1) {
	        run("Analyze Angles", "buffer=0 window=250 overlap=0.25 save_path='C:\\Users\\gobes001\\LocalSoftware\\AnalyseDirectionality\\Python Scripts\\input\\testing' python_path='py' script_path='C:\\Users\\gobes001\\LocalSoftware\\AnalyseDirectionality\\Python Scripts' run_python=True python_arguments='-r'");
			close();
		}
}
setBatchMode(false);
