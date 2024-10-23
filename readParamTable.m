% function [cellData,rev] = readParamTable(filename,sheet)  
%
%     Loads Excel spreadsheet containing parameters that define a specific
%     lithium-ion cell, and stores these parameters in a structure form.
%
%   Inputs:
%     filename: Filename for Excel spreadsheet holding parameter values
%     sheet: Name of worksheet in the spreadsheet that should be loaded
%
%   Outputs:     
%     cellData: A structure containing the params defined in spreadsheet
%     rev:      The revision number of the spreadsheet format.

% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.

function [cellData,rev] = readParamTable(fileName,sheet)
  if ~exist(fileName,'file'),
    error('Spreadsheet "%s" does not exist.',fileName);
  end
  [~,sheets] = xlsfinfo(fileName);
  if sum(strcmpi(sheets,sheet)) == 0,
    error('Spreadsheet "%s" does not contain worksheet "%s".',fileName,sheet);
  end  
  [~,~,data]  = xlsread(fileName,sheet); % Read Excel file

  % Extract version from data
  rev = data{3,11};
  
  % Scan file for markers delineating certain sections.
  % These sections should occur in this order in the Excel spreadsheet.
  for ii = 1:length(data)
    switch lower( data{ii,1} )
      case {'cell information'},         yName = ii;
      case {'constants','constant'},     yConst = {ii; 'const'};
      case {'negative electrode','neg'}, yNeg = {ii; 'neg'};
      case {'seperator','sep'},          ySep = {ii;'sep'};
      case {'positive electrode','pos'}, yPos = {ii;'pos'};
      otherwise, % do nothing
    end
  end

  ind = [yConst,yNeg,ySep,yPos,{length(data);0}];
  % Strip off all NaN cells
  for ii = 1:(length(ind)-1)
    for jj = (ind{1,ii}+2):(ind{1,ii+1}-1)
      param{jj,ii} = data{jj,2}; %#ok<AGROW>
      value{jj,ii} = data{jj,3}; %#ok<AGROW>
      % unit{jj,ii}  = data{jj,4};  Not used now, but maybe later.
      if( isnan( param{jj,ii} ) ) 
        param{jj,ii} = [];  %#ok<AGROW>
      end
      if( isnan( value{jj,ii} ) ) 
        value{jj,ii} = []; %#ok<AGROW>
      end
      % if( isnan( unit{jj,ii}  ) ) unit{jj,ii}  = []; end
    end
  end

  % Form the region structs
  name = data{yName+2,3};
  for ii = 1:(length(ind)-1)
    eval(sprintf('%s = struct();',ind{2,ii}));
    eval(sprintf('%s = setfield(%s,''name'',''%s'');',...
                 ind{2,ii},ind{2,ii},ind{2,ii}));
    for jj = 1:length(param)
     if(~isempty( param{jj,ii} ))
       eval(sprintf('%s = setfield(%s,param{jj,ii},value{jj,ii});',...
                 ind{2,ii},ind{2,ii}));
     end
    end
  end
  
  % Put everything together into cellData structure
  cellData = struct('name',name,'const',const,'neg',neg,'sep',sep,'pos',pos);
end