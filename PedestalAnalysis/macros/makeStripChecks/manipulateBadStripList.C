void manipulateBadStripList (string inputListTXT){

  TString outputFile (inputListTXT.c_str());
  outputFile.ReplaceAll(".txt","_modified.txt");
  ofstream output(outputFile.Data());

  long int badStrips = 0;

  ifstream file (inputListTXT.c_str());
  if(file.is_open()){
    string line;
    while(!file.eof()){
      getline(file,line);
      stringstream ss (line);
      vector<string> column;
      if(line == "") continue;
      while(ss >> line)
      	column.push_back(line);
      if(column.size() != 2)
      	cerr<<"Problem with the text file input "<<endl;
      if(column.back() == "0")
      	continue;
      else{
      	output << ss.str() <<"\n";
	badStrips = badStrips + atof(column.back().c_str());
      }
    }    
  }

  cout<<"Total number bad strips "<<badStrips<<endl;
  cout<<"Output File "<<outputFile<<endl;

  file.close();
  output.close();

}
