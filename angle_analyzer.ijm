input = "D:\\Data\\Microscopy\\Airyscan\\2023\\DP1\\Test\\";
setBatchMode(true);
list = getFileList(input);
for (i = 0; i < list.length; i++){
		open(input + list[i]);
		if (nImages>=1) {
	        run("Analyze Angles", "buffer=0 window=350 overlap=0.25 save_path='C:\\Users\\gobes001\\LocalSoftware\\AnalyseDirectionality\\Python Scripts\\input\\testing' python_path='py' script_path='C:\\Users\\gobes001\\LocalSoftware\\AnalyseDirectionality\\Python Scripts' run_python=True");
			close();
		}
}
setBatchMode(false);
