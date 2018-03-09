# PeptoPred

2.16.18 Day 1:
Learned to use git on shell. Is it wise to commit new project file structure to git? For now, just committing diary to git repo.

2.18.18 Day 3:
Learned to do some one-liners and use wildcards (star, questionmark) on shell. Maybe I can use them to commit all files in repo at once.

2.19.18 Day 4:
Starting the switch to Sublime Text! Was able to enter link to Sublime.app in PATH for easy terminal access. Also created little folder maker script, but for what purpose? After reading the quick guide to computational bio, I have many more questions than answers. What is a make file? Maybe this folder maker script will be useful once we have makefiles, specific calls to parsers and more.

3.2.18 Day ?:
Created a very ugly encoder using a hard-coded feature map and a few bad if-for statements. I wonder if we could use binary to encode each amino acid using only 5 1/0's - would only take 5 positions instead of 20. First it's more important to test feeding into svm or randomforest. Enjoyed the whale talk today. My favorite blog is http://infoproc.blogspot.com/

3.6.18 Day ??:
Much work was done this day. Successfully ran my first prediction through a model trained on ~10 proteins - testpeps.txt, of which 2 or 3 had signal peptides. The model missed both signal peptides in the test dataset - testpeps2.txt, just classified everything as globular. Haven't written output yet - just need to save a list of the proteins and their encodings as they are encoded - should be able to write something up to call them after prediction for pred formatting. Read today about featuretools.com. First want to try PSSM. Reading on SignalP's implementation - downloaded their desktop version of the tool but it's linux only - maybe will try looking at their code in the computer lab. Also tried creating a model based on the full dataset but aborted after running for an hour. Started thinking about program organization and re-organized formatting.py - now global sliding window gen is in "FormatSansDict.py"

3.9.18 Day ???:
Full model generated from LinearSVC. Plan to run several window sizes and algorithm types overnight. First need to test if joblib.dump/load 'wb' works for features.
