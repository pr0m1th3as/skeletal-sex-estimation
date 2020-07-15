% Copyright (C) 2020 Andreas Bertsatos <abertsatos@biol.uoa.gr>
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, see <http://www.gnu.org/licenses/>.

% A function for skeletal sex estimation.
%
% The function estimate_sex provides a graphical user interface for estimating
% sex from various sleletal elements according to the selected classifier saved
% in the respective .mat file. The classifier files can be parsed as an input
% argument to the function by invoking
%
% e.g. estimate_sex("long_bone_classifier.mat")
%
% or can be selected through the open file dialog when the function is called
% withough any input arguments.
%
% Currently, the function supports two types of classifiers.
%
% Use "long_bone_classifier.mat" for estimating sex from CSG data of the femur,
% tibia, and humerus bones. The CSG data, which must conform to the same format
% as in the "CSGtestdata.csv" file, can be automatically extracted from virtual
% 3D bone models with the long-bone-diaphyseal-CSG-Toolkit, which also saves
% the required data to the appropriate format. The "long_bone_classifier.mat"
% bundled with this function utilizes statistical models derived from a modern
% Greek population skeletal reference collection, the Athens Collection. 
%
% Use "vertebrae_classifier.mat" for estimating sex from linear measurements of
% the T1, T12 and L1 vertebrae. This classifier utilizes discriminant functions
% derived from two European populations, a modern Greek (Athens Collection)
% and a contermporary Danish (University of Copenhagen CT-scan data).
function estimate_sex(varargin)
  % check if 'io' package is loaded (load it if necessary)
  if ~exist("csv2cell")
    pkg load io
  endif
  % initial check for input arguments
  if nargin == 0
    class_fname = uigetfile({"*.mat", "data text container"},...
                             "Select the appropriate classifier container file");
  elseif nargin == 1
    class_fname = varargin{1};
  else
    error("Wrong number of input arguments\n");
  endif
  % get the data target group of the classifier from the 'type' filed of the 'description' structure
  classifier = load(class_fname);
  if ~isfield(classifier, "classifiers") || ~isfield(classifier, "description")
    printf("The selected file is not a valid classifier container!\n"); return;
  endif
  type = classifier.description.datatype;
  % check the type of classifier to request for specific type of data set
  tt = find(strcmp({type}, {"CSG-Toolkit", "skullanalyzer", "Vertebrae"}));
  switch (tt)
    % for CSG-Toolkit data
    case 1
      tdata_fname = uigetfile({"*.csv", "longbone CSG properties"},...
                              "Select a dataset created with the CSG-Toolkit");
      % read dataset from csv file and check its integrity
      all_data = csv2cell(tdata_fname);
      if (size(all_data,1) < 1 || size(all_data,2) ~= 47)
        printf("error: The selected file is not a valid CSG properties container!\n\
       Please, use the inspectCSG.m function from the CSG-Toolkit to\n\
       create an appropriately formated csv file with the CSG properties\n\
       of the long bones, which you wish to estimate sex for.\n"); return;
      endif
      % remove the header row and keep a list of sample names from the first column
      sample_IDs = all_data([2:end],1);
      clear all_data;
      % reload the csv file in a matrix, keep only the necessary CSG data, calculate
      % the extra variables used for sex estimation and rearrange the matrix accordingly
      CSG = csvread(tdata_fname);
      CSG(1,:) = []; CSG(:,[1,3:7,12,15,20,23,28,31,36,39,44,47]) = [];
      % keep maxDistance, Areas, Perimeters, Ix, Iy, Imin, Imax and create new variables with ratios
      % of Ix/Iy and Imax/Imin as well as ArPerIndex for each cross section
      CSGnV = CSG(:,[1:3]);
      CSGnV = [CSGnV, ((CSG(:,2))*(4*pi))./(CSG(:,3).^2), CSG(:,[4:5]), CSG(:,4)./CSG(:,5)];
      CSGnV = [CSGnV, CSG(:,[6:7]), CSG(:,7)./CSG(:,6), CSG(:,[8:9])];
      CSGnV = [CSGnV, ((CSG(:,8))*(4*pi))./(CSG(:,9).^2), CSG(:,[10:11]), CSG(:,10)./CSG(:,11)];
      CSGnV = [CSGnV, CSG(:,[12:13]), CSG(:,13)./CSG(:,12), CSG(:,[14:15])];
      CSGnV = [CSGnV, ((CSG(:,14))*(4*pi))./(CSG(:,15).^2), CSG(:,[16:17]), CSG(:,16)./CSG(:,17)];
      CSGnV = [CSGnV, CSG(:,[18:19]), CSG(:,19)./CSG(:,18), CSG(:,[20:21])];
      CSGnV = [CSGnV, ((CSG(:,20))*(4*pi))./(CSG(:,21).^2), CSG(:,[22:23]), CSG(:,22)./CSG(:,23)];
      CSGnV = [CSGnV, CSG(:,[24:25]), CSG(:,25)./CSG(:,24), CSG(:,[26:27])];
      CSGnV = [CSGnV, ((CSG(:,26))*(4*pi))./(CSG(:,27).^2), CSG(:,[28:29]), CSG(:,28)./CSG(:,29)];
      CSGnV = [CSGnV, CSG(:,[30:31]), CSG(:,31)./CSG(:,30)];
      
      % call the main ui for sex estimation based on CSG data
      h = user_interface_CSG(classifier, sample_IDs, CSGnV);
      
    % for skullanalyzer data      
    case 2
      tdata_fname = uigetfile({"*.mat", "skullanalyzer morphometric features"},...
           "Select one or multiple .mat files created with the skullanalyzer", "MultiSelect", "on");
      % check that all selected .mat files contain skullanalyzer's morphometric features,
      % dismiss those that don't, and create a list of sample names with the valid filenames.
      % by removing the .mat extension
      valid_samples = 0;
      % if 1 file selected, convert char string to a cell string
      if ischar(tdata_fname)
        tdata_fname = {tdata_fname};
      else
        tdata_fname = tdata_fname';
      endif
      for i=1:length(tdata_fname)
        vname = tdata_fname{i};
        if strcmp(".mat", vname([end-3:end]))
          temp = load(tdata_fname{i});
          if isfield(temp, "EFD") && isfield(temp, "HMI") && isfield(temp, "FCC")
            valid_samples += 1;
            sample_IDs(valid_samples) = {sprintf("%s", vname([1:end-4]))};
          end
          clear temp;
        endif
      endfor
      
      if valid_samples == 0
        printf("error: The selected file(s) is not a valid morphometric features container!\n\
       Please, use the skullanalyzer to create appropriately formated .mat file(s)\n\
       from your cranial sample(s), which you wish to estimate sex for.\n"); return;
      elseif valid_samples < i
        printf("Some files are not valid morphometric features container and were excluded.\n");
      end
      % transpose sample list 
      sample_IDs = sample_IDs';
      
      % call the main ui for sex estimation based on CSG data
      h = user_interface_skullanalyzer(classifier, sample_IDs);
      
    % for vertebral linear measurement data      
    case 3
      % call the main ui for sex estimation based on vertebral linear measurements
      h = user_interface_vertebrae(classifier);
      
    otherwise
      printf("error: The selected file is not a valid classifier container!\n"); return;
  endswitch
endfunction

%
% main user interface for sex estimation based on CSG data
function h = user_interface_CSG(classifier, sample_IDs, CSGdata)
  % list bones and sides for the CSG-Toolkit sex classifier
  bone = {"Femur", "Tibia" "Humerus"};
  side = {"Left", "Right"};
  % UI background colour
  bg_color = [0.8 0.8 0.8];
  % create a main window for user interface and populate it with necessary user options
  h.f = figure("Name", "skeletal sex estimation based on CSG-Toolkit data", "numbertitle", "off",...
               "menubar", "none", "color", bg_color, "Position", [300, 300, 560, 450],...
               "resize", "off");
  % display available classification methods (RBF or LDA)
  class_methods = {"RBF kernel SVM classification"; "Linear Discriminat Function Analysis"};
  h.class_method = uibuttongroup("title", "Available classification methods", "titleposition",...
                                 "centertop", "backgroundcolor", bg_color,...
                                 "position", [0.2, 0.76, 0.6, 0.2]);
    uicontrol(h.class_method, "style", "radiobutton", "string", class_methods{1}, "tag", "RBF",...
                              "backgroundcolor", bg_color, "position", [50, 40, 300, 30]);
    uicontrol(h.class_method, "style", "radiobutton", "string", class_methods{2}, "tag", "LDA",...
                              "backgroundcolor", bg_color, "position", [50, 10, 300, 30]);
  set(get(h.class_method, "children")(2), "selected", "on", "value", 1);
  % selection of skeletal sample, bone and side (if applicable)
  h.text_S = uicontrol(h.f, "style", "text", "string", "Available samples",...
                            "backgroundcolor", bg_color, "position", [30, 300, 200, 30]);
  h.sample = uicontrol(h.f, "style", "popupmenu", "string", sample_IDs,...
                            "position", [30, 270, 200, 30]);
  h.text_B = uicontrol(h.f, "style", "text", "string", "Select bone",...
                            "backgroundcolor", bg_color, "position", [280, 300, 100, 30]);
  h.bone_type = uicontrol(h.f, "style", "popupmenu", "string", bone,...
                               "position", [280,270,100,30]);
  h.text_S = uicontrol(h.f, "style", "text", "string", "Select side",
                            "backgroundcolor", bg_color, "position", [430, 300, 100, 30]);
  h.bone_side = uicontrol(h.f, "style", "popupmenu", "string", side,...
                               "position", [430, 270, 100, 30]);
  % selection of particular classifiers
  h.class1 = uicontrol(h.f, "style", "checkbox", "string", "Classifier #1:", "fontangle",...
                            "italic", "fontweight", "bold", "backgroundcolor", bg_color,...
                            "position", [30, 210, 120, 20]);
  h.class2 = uicontrol(h.f, "style", "checkbox", "string", "Classifier #2:", "fontangle",...
                            "italic", "fontweight", "bold", "backgroundcolor", bg_color,...
                            "position", [30, 160, 120, 20]);
  h.class3 = uicontrol(h.f, "style", "checkbox", "string", "Classifier #3:", "fontangle",...
                            "italic", "fontweight", "bold", "backgroundcolor", bg_color,...
                            "position", [30, 110, 120, 20]);
  set(h.class1, "value", 1, "selected", "on");
  % containers for each classifier's results
  results_string = {"","",""};
  h.result1 = uicontrol(h.f, "style", "text", "string", results_string{1}, "horizontalalignment",...
                             "left", "backgroundcolor", bg_color, "position", [160, 210, 400, 20]);
  h.result2 = uicontrol(h.f, "style", "text", "string", results_string{2}, "horizontalalignment",...
                             "left", "backgroundcolor", bg_color, "position", [160, 160, 400, 20]);
  h.result3 = uicontrol(h.f, "style", "text", "string", results_string{3}, "horizontalalignment",...
                             "left", "backgroundcolor", bg_color, "position", [160, 110, 400, 20]);
  % push buttons for estimating sex and saving the results in a csv file
  estimate = uicontrol (h.f, "style", "pushbutton", "string", "estimate sex",...
                             "Position", [130, 20, 100, 50], "backgroundcolor", bg_color);
  savedata = uicontrol (h.f, "style", "pushbutton", "string", "save results",...
                             "Position", [340, 20, 100, 50], "backgroundcolor", bg_color);
  % set callbacks
  set (estimate, "callback", {@run_estimate_CSG, estimate, h, classifier, CSGdata});
  set (savedata, "callback", {@run_savedata_CSG, savedata, h, classifier, CSGdata});
endfunction

% callback function for estimating sex from CSG data according to user's selection
function run_estimate_CSG(estimate, dummy1, dummy2, h, classifier, CSGdata)
  % take current user selections and create indexing accordingly
  [method, sample_idx, sample_ID, bone, side, skelement, DFs] = update_CSG_UI(h);
  % get sample from CSGdata
  sample_data = CSGdata(sample_idx,:);
  % for each DF keep the appropriate variables, normalize them and estimate sex
  % and posterior probability according to the selected method
  [sample_sex, sample_pb] = predict_sex_CSG(DFs, method, skelement, classifier, sample_data);
  % create the results' strings for each selected DF and update the ui
  r_idx = 0;
  for c_idx = DFs
    r_idx += 1;
    result_str = sprintf("Sample is %s with posterior probability %0.2f",...
                          sample_sex{r_idx}, sample_pb(r_idx));
    r_sel = cstrcat("result",sprintf("%i",c_idx));
    set(h.(r_sel), "string", result_str);
  endfor
endfunction

% callback function for saving results from CSG data according to user's selection
function run_savedata_CSG(savedata, dummy1, dummy2, h, classifier, CSGdata)
  % take current user selections and create indexing accordingly
  [method, sample_idx, sample_ID, bone, side, skelement, DFs] = update_CSG_UI(h);
  % get sample from CSGdata
  sample_data = CSGdata(sample_idx,:);
  % for each DF keep the appropriate variables, normalize them and estimate sex
  % and posterior probability according to the selected method
  [sample_sex, sample_pb] = predict_sex_CSG(DFs, method, skelement, classifier, sample_data);
  % create a header and a standard filename for the results' file
  results_header = {"Sample ID", "bone", "side", "method", "classifier",...
                    "predicted sex", "posterior probability"};
  rfname = "estimate_sex results.csv";
  % check if file already exists (if not, create one with header only)
  if (~exist(rfname))
    cell2csv(rfname, results_header);
  endif
  % read file and check its header size (in case it existed already)
  results = csv2cell(rfname);
  if (size(results, 2) ~= 7)
    results = results_header;
  endif
  % append a new entries at the end row with current sex estimation results
  r_idx = 0;
  for c_idx = DFs
    r_idx += 1;
    results(end+1,:)  = {sample_ID, bone, side, method, c_idx, sample_sex{r_idx}, sample_pb(r_idx)};
  endfor
  cell2csv(rfname, results);
endfunction

% for each selected DF keep the appropriate variables, normalize them and estimate sex
% and posterior probability according to the selected method regarding CSG data
function [sample_sex, sample_pb] = predict_sex_CSG(DFs, method, skelement, classifier, sample_data)
  if length(DFs) > 0
    for i=1:length(DFs)
      % get index of classifier from description table
      temp_dt = classifier.description.(method)();
      classifier_idx = cell2mat(temp_dt(find(strcmp(skelement, temp_dt(:,1))),DFs(i)+1));
      clear temp_dt;
      % get normalization coefficients and variable indices for corresponding selected classifier
      variable_indices = classifier.normalcoef(classifier_idx).var_num;
      norm_mean = classifier.normalcoef(classifier_idx).mean;
      norm_StD = classifier.normalcoef(classifier_idx).StD;
      % truncate to required variables and normalize data
      sample_trunc_data = sample_data(variable_indices);
      sample_norm_data = (sample_trunc_data - norm_mean) ./ norm_StD;
      % get appropriate classifier
      sel_class = classifier.classifiers(classifier_idx).(method);
      % evaluate sample's data with classifier
      df(i) = evaluate_model(sample_norm_data, sel_class);
      % some cleanup
      clear variable_indices sample_norm_data sample_trunc_data norm_mean norm_StD sel_class
      % evaluate sample's posterior probability of the calculated outcome
      temp_pb = classifier.post_proba(classifier_idx).(method);
      if isfield(temp_pb, "sectioning_point") % for LDA classifier
        sp = temp_pb.sectioning_point;
        cf = temp_pb.centroid_females;
        cm = temp_pb.centroid_males;
        % find distance from sectioning point
        dist = df(i) - sp;
        if cf < cm    % female centroid smaller than male centroid
          if dist < 0 % then female
            sample_sex(i) = {"female"};
            % calculate posterior probability for sample belonging to female group
            sample_pb(i) = exp(-dist)/(exp(-dist)+exp(dist));
           else       % then male
            sample_sex(i) = {"male"};
            % calculate posterior probability for sample belonging to female group
            sample_pb(i) = exp(dist)/(exp(dist)+exp(-dist));
          endif
        else          % male centroid smaller than female centroid
          if dist < 0 % then male
            sample_sex(i) = {"male"};
            % calculate posterior probability for sample belonging to female group
            sample_pb(i) = exp(-dist)/(exp(-dist)+exp(dist));
           else       % then female
            sample_sex(i) = {"female"};
            % calculate posterior probability for sample belonging to female group
            sample_pb(i) = exp(dist)/(exp(dist)+exp(-dist));
          endif
        endif
        % some cleanup
        clear sp cf cm
      elseif isfield(temp_pb, "discrete_PDF") % for RBF classifier
        dPDF = temp_pb.discrete_PDF;
        fg = temp_pb.female_group;
        mg = temp_pb.male_group;
        % get posterior probability according to the ranges provided in the lookup table
        sample_pb(i) = dPDF(find(dPDF(:,1) <= abs(df(i)) & dPDF(:,2) > abs(df(i))),3);
        % check if first group (negative RBF scores) is males or females
        if mg < fg    % male group takes negative scores (females have positive values)
          if df(i) < 0 % then male
            sample_sex(i) = {"male"};
           else       % then female
            sample_sex(i) = {"female"};
          endif
        else          % female group takes negative scores (males have positive values)
          if df(i) < 0 % then female
            sample_sex(i) = {"female"};
           else       % then male
            sample_sex(i) = {"male"};
          endif
        endif
        % some cleanup
        clear dPDF fg mg
      endif
      % some cleanup
       clear temp_pb
    endfor
  endif
endfunction

% take current user selections regarding CSG data and create indexing accordingly
function [method, sample_idx, sample_ID, bone, side, skelement, DFs] = update_CSG_UI(h)
  % get classification method
  method = get(get(h.class_method, "selectedobject"), "tag");
  % get selected sample
  sample_str = get(h.sample, "string");
  sample_idx = get(h.sample, "value");
  sample_ID = sample_str{sample_idx};
  % get selected bone
  bone_str = get(h.bone_type, "string");
  bone_idx = get(h.bone_type, "value");
  bone = bone_str{bone_idx};
  % get selected side
  side_str = get(h.bone_side, "string");
  side_idx = get(h.bone_side, "value");
  side = side_str{side_idx};
  % concatenate bone and side for searching the index table of the classifier's description
  skelement = cstrcat(bone, " ", side);
  % get selected classifiers
  DF1 = get(h.class1, "value");
  DF2 = get(h.class2, "value");
  DF3 = get(h.class3, "value");
  DFs = [1, 2, 3];
  DFs = DFs(find(DFs .* [DF1, DF2, DF3]));
endfunction



%
% main user interface for sex estimation based on vertebral linear measurements
function h = user_interface_vertebrae(classifier)
  % get description data from classifier
  vertebra = classifier.description.vertebra;
  population = classifier.description.population;
  % UI background colour
  bg_color = [0.8 0.8 0.8];
  % create a main window for user interface and populate it with necessary user options
  h.f = figure ("Name", "skeletal sex estimation based on vertebral linear measurements",...
                "numbertitle", "off", "menubar", "none", "color", bg_color, "Position",...
                [300 300 560 450], "resize", "off");
  % display available classification methods (LDA only)
  class_methods = {"Linear Discriminat Function Analysis"};
  h.class_method = uibuttongroup("title", "Available classification methods", "titleposition",...
                                 "centertop", "backgroundcolor", bg_color,...
                                 "position", [0.2, 0.82, 0.6, 0.14]);
    uicontrol(h.class_method, "style", "radiobutton", "string", class_methods{1}, "tag", "LDA",...
                              "backgroundcolor", bg_color, "position", [40, 10, 300, 30]);
  set(get(h.class_method, "children")(1), "selected", "on", "value", 1);
  % selection of vertebra and population
  h.text_V = uicontrol(h.f, "style", "text", "string", "Select vertebra",...
                            "backgroundcolor", bg_color, "position", [30, 325, 110, 30]);
  h.vertebra = uicontrol(h.f, "style", "popupmenu", "string", vertebra,...
                              "position", [30, 295, 110, 30]);
  h.text_P = uicontrol(h.f, "style", "text", "string", "Select population",...
                            "backgroundcolor", bg_color, "position", [180, 325, 150, 30]);
  h.population = uicontrol(h.f, "style", "popupmenu", "string", population,...
                                "position", [180, 295, 150, 30]);
  % text container for sample ID
  h.text_S = uicontrol(h.f, "style", "text", "string", "Sample ID",...
                            "backgroundcolor", bg_color, "position", [380, 325, 150, 30]);
  h.sample = uicontrol(h.f, "style", "edit", "string", "", "position", [380, 295, 150, 30]);
  % display measurement input boxes for the default classifier
  % find default classifier and corresponding number of measurements 
  skelement = cstrcat(population{1}, " ", vertebra{1});
  idx = find(strcmp(skelement, classifier.description.LDA(:,1)));
  class_idx = classifier.description.LDA(idx, 2){};
  measurements = classifier.classifiers(class_idx).LDA.measurements;
  varnum = length(measurements);
  for i = 1:4
    if i <= varnum
      h.var(i).name = uicontrol(h.f, "style", "text", "string", measurements{i}, "fontangle",...
                                     "italic", "horizontalalignment", "left", "backgroundcolor",...
                                     bg_color, "position", [30, 255 - ((i - 1) * 40), 260, 20]);
      h.var(i).value = uicontrol(h.f, "style", "edit", "string", "",...
                                               "position", [300, 255 - ((i - 1) * 40), 50, 20]);
    else
      h.var(i).name = uicontrol(h.f, "style", "text", "string", "", "fontangle", "italic",...
                                     "horizontalalignment", "left", "backgroundcolor", bg_color,...
                                     "position", [30, 255 - ((i - 1) * 40), 260, 20]);
      h.var(i).value = uicontrol(h.f, "style", "edit", "string", "", "visible", "off",...
                                      "position", [300, 255 - ((i - 1) * 40), 50, 20]);
    endif
  endfor
  % prepare a container for the classifier's results
  results_string = {""};
  h.results = uicontrol(h.f, "style", "text", "string", results_string{1}, "horizontalalignment",...
                             "center", "backgroundcolor", bg_color, "position", [80, 100, 420, 20]);
  % push button for displaying measuremenents, estimating sex, and saving the results in a csv file
  enterval = uicontrol (h.f, "style", "pushbutton", "string", "display\n measurements",...
                             "Position", [70, 20, 120, 50], "backgroundcolor", bg_color);
  estimate = uicontrol (h.f, "style", "pushbutton", "string", "estimate sex",...
                             "Position", [240, 20, 100, 50], "backgroundcolor", bg_color);
  savedata = uicontrol (h.f, "style", "pushbutton", "string", "save results",...
                             "Position", [390, 20, 100, 50], "backgroundcolor", bg_color);
  % set callbacks
  set (enterval, "callback", {@run_enterval_vertebrae, enterval, h, classifier});
  set (estimate, "callback", {@run_estimate_vertebrae, estimate, h, classifier});
  set (savedata, "callback", {@run_savedata_vertebrae, savedata, h, classifier});
endfunction

% callback function for updating the vertebral measurements according to user's selection
function run_enterval_vertebrae(enterval, dummy1, dummy2, h, classifier)
  % take current user selections and display the corresponding measurements for user input
  % get selected vertebra
  vertebra_str = get(h.vertebra, "string");
  vertebra_idx = get(h.vertebra, "value");
  vertebra = vertebra_str{vertebra_idx};
  % get selected population
  population_str = get(h.population, "string");
  population_idx = get(h.population, "value");
  population = population_str{population_idx};
  % concatenate bone and side for searching the index table of the classifier's description
  skelement = cstrcat(population, " ", vertebra);
  % get selected classifier and display the corresponding measurements and input boxes
  idx = find(strcmp(skelement, classifier.description.LDA(:,1)));
  class_idx = classifier.description.LDA(idx, 2){};
  measurements = classifier.classifiers(class_idx).LDA.measurements;
  varnum = length(measurements);
  for i = 1:4
    if i <= varnum
      set(h.var(i).name, "string",  measurements{i});
      set(h.var(i).value, "visible", "on");
    else
      set(h.var(i).name, "string",  "");
      set(h.var(i).value, "visible", "off");
    endif
  endfor
endfunction

% callback function for estimating sex from user defined vertebral measurements
function run_estimate_vertebrae(estimate, dummy1, dummy2, h, classifier)
  % get sample's ID, sex, posterior probability, used measurements and their values
  [sample_ID, sample_sex, sample_pb, M_names, S_data] = predict_sex_vertebrae(h, classifier);
  % create a string for displaying classification results and update the ui
  result_str = sprintf("Sample %s is %s with posterior probability %0.2f",...
                       sample_ID, sample_sex{:}, sample_pb);
  set(h.results, "string", result_str);
endfunction

% callback function for saving sex estimation results from vertebral measurements
function run_savedata_vertebrae(savedata, dummy1, dummy2, h, classifier)
  % get sample's ID, sex, posterior probability, used measurements and their values
  [sample_ID, sample_sex, sample_pb, M_names, S_data] = predict_sex_vertebrae(h, classifier);
  % pad measurement name and sample data vectors with blanks
  varnum = length(M_names);
  S_data = num2cell(S_data);
  for i = 1:4
    if i > varnum
      M_names(i) = {""};
      S_data(i) = {""};
    endif
  endfor
  % get selected vertebra
  vertebra_str = get(h.vertebra, "string");
  vertebra_idx = get(h.vertebra, "value");
  vertebra = vertebra_str{vertebra_idx};
  % get selected population
  population_str = get(h.population, "string");
  population_idx = get(h.population, "value");
  population = population_str{population_idx};
  % create a header and a standard filename for the results' file
  results_header = {"Sample ID", "population", "vertebra", "predicted sex",...
                    "posterior probability", "measurement", "value", "measurement", "value",...
                    "measurement", "value", "measurement", "value"};
  rfname = "estimate_sex results.csv";
  % check if file already exists (if not, create one with header only)
  if (~exist(rfname))
    cell2csv(rfname, results_header);
  endif
  % read file and check its header size (in case it existed already)
  results = csv2cell(rfname);
  if (size(results, 2) ~= 13)
    results = results_header;
  endif
  % append a new entry at the end row with current sex estimation results
  results(end+1,:)  = {sample_ID, population, vertebra, sample_sex{:}, sample_pb,...
                       M_names{1}, S_data{1}, M_names{2}, S_data{2}, M_names{3}, S_data{3},...
                       M_names{4}, S_data{4}};
  % save to csv file
  cell2csv(rfname, results);
endfunction

% get selected options and classify a sample from vertebral measurements
% return sample's ID, sex, posterior probability, used measurements and their values
function [sample_ID, sample_sex, sample_pb, M_names, S_data] = predict_sex_vertebrae(h, classifier)
  % get selected vertebra
  vertebra_str = get(h.vertebra, "string");
  vertebra_idx = get(h.vertebra, "value");
  vertebra = vertebra_str{vertebra_idx};
  % get selected population
  population_str = get(h.population, "string");
  population_idx = get(h.population, "value");
  population = population_str{population_idx};
  % concatenate bone and side for searching the index table of the classifier's description
  skelement = cstrcat(population, " ", vertebra);
  % find corresponding classifier and copy DF weights into neew structure
  idx = find(strcmp(skelement, classifier.description.LDA(:,1)));
  class_idx = classifier.description.LDA(idx, 2){};
  DF.weights = classifier.classifiers(class_idx).LDA.weights;
  % retrieve sample's measurements from UI
  M_names = classifier.classifiers(class_idx).LDA.measurements;
  varnum = length(M_names);
  S_data = [];
  for i = 1:4
    if i <= varnum
      S_data = [S_data, str2num(get(h.var(i).value, "string"))];
    endif
  endfor
  % evaluate sample's data with corresponding classifier
  df = evaluate_model(S_data, DF);
  % get sectioning point and group centroids to classify sample and compute posterior probability
  sp = classifier.post_proba(class_idx).LDA.sectioning_point;
  cf = classifier.post_proba(class_idx).LDA.centroid_females;
  cm = classifier.post_proba(class_idx).LDA.centroid_males;
  % find distance from sectioning point
  dist = df - sp;
  if cf < cm    % female centroid smaller than male centroid
    if dist < 0 % then female
      sample_sex = {"female"};
      % calculate posterior probability for sample belonging to female group
      sample_pb = exp(-dist)/(exp(-dist)+exp(dist));
     else       % then male
      sample_sex = {"male"};
      % calculate posterior probability for sample belonging to female group
      sample_pb = exp(dist)/(exp(dist)+exp(-dist));
    endif
  else          % male centroid smaller than female centroid
    if dist < 0 % then male
      sample_sex = {"male"};
      % calculate posterior probability for sample belonging to female group
      sample_pb = exp(-dist)/(exp(-dist)+exp(dist));
     else       % then female
      sample_sex = {"female"};
      % calculate posterior probability for sample belonging to female group
      sample_pb = exp(dist)/(exp(dist)+exp(-dist));
    endif
  endif
  % get user defined sample ID
  sample_ID = get(h.sample, "string");
endfunction


% function to evaluate LDA and RBF kernel based SVM models
function df = evaluate_model(data, classifier)
	% determine the type of classifier by the number of fields
	if(length(fieldnames(classifier))==4)
		% apply RBF
		df = 0;
		for i=1:length(classifier.dual_coef)
			d(i) = norm(classifier.support_vector(i,:) - data);
			t(i) = exp(-classifier.gamma_param * (d(i) ^ 2));
			df += classifier.dual_coef(i) * t(i);
		endfor
		df += classifier.rho_intercept;		
	elseif(length(fieldnames(classifier))==1)
		% apply LDA (constant + coefficients x measurements)
		df = classifier.weights(1) + sum(classifier.weights([2:end]) .* data);
	endif
endfunction