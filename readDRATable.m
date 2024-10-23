% function [draData,tflist,rev] = readDRATable(fileName,sheet)  
%
%     Loads Excel spreadsheet containing DRA settings, and converts the
%     settings to structure form.
%
%   Inputs:
%     fileName: Filename for Excel spreadsheet holding DRA settings
%     sheet:    Name of worksheet in the spreadsheet that should be loaded
%
%   Outputs:     
%     draData:  A structure containing the params defined in spreadsheet
%     tflist:   List of transfer functions to be evaluated by DRA
%     rev:      The revision number of the spreadsheet format.

% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.

function [draData,tflist,rev] = readDRATable(fileName,sheet) 
  if ~exist(fileName,'file'),
    error('Spreadsheet "%s" does not exist.',fileName);
  end
  [~,sheets] = xlsfinfo(fileName);
  if sum(strcmpi(sheets,sheet)) == 0,
    error('Spreadsheet "%s" does not contain worksheet "%s".',fileName,sheet);
  end  
  [~,~,data]  = xlsread(fileName,sheet); % Read Excel file

  % Extract version from data
  rev = data{3,9};
  
  % Scan file for markers delineating certain sections.
  % These sections should occur in this order in the Excel spreadsheet.
  for ii = 1:length(data),
    switch lower( data{ii,1} )
      case {'cell information'},               yName = ii; %#ok<NASGU>
      case {'environmental'},                  yEnv = {ii; 'env'};
      case {'dra'},                            yDRA = {ii; 'dra'};
      case {'transfer function','tf','tf''s'}, yTF  = {ii; 'tf'};
      otherwise, % do nothing
    end
  end

  ind = [yEnv,yDRA,yTF,{length(data);0}];
  % Strip Off All NaN Cells
  for ii = 1:(length(ind)-1),
    for jj = (ind{1,ii}+2):(ind{1,ii+1}-1),
      param{jj,ii} = data{jj,2}; %#ok<AGROW>
      value{jj,ii} = data{jj,3}; %#ok<AGROW>
      % unit{jj,ii}  = data{jj,4}; Not used right now, but maybe later
      if( isnan( param{jj,ii} ) ) 
        param{jj,ii} = [];  %#ok<AGROW>
      end
      if( isnan( value{jj,ii} ) ) 
        value{jj,ii} = [];  %#ok<AGROW>
      end
      % if( isnan( unit{jj,ii}  ) ) unit{jj,ii}  = []; end
    end
  end

  % Form the DRA Parameters
  for ii = 1:(length(ind)-1)
    for jj = 1:length(param)
      % For everything but the transfer functions
      if(~isempty( param{jj,ii} ) && (jj < yTF{1}) )
        if( ischar( value{jj,ii} ) )
          num = eval(value{jj,ii});
        else
          num = value{jj,ii};
        end
        draData.(param{jj,ii}) = num;
      end

      % For transfer function list
      if( strcmpi(param{jj,ii},'tflist') )
        kk = 1;
        tflist = {sprintf('%s',value{jj,ii})};
        while ( ~isempty( value{jj+kk,ii} ) ),
          tflist = [tflist; {sprintf('%s',value{jj+kk,ii})} ]; %#ok<AGROW>
          kk = kk+1;
        end
      end
    end
  end
end