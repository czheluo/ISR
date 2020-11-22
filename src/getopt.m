function [opt unprocessed_arg] = getopt(opt,varargin)
% To be used in functions for parsing varargin (the list of input arguments of a variable length)
% and extracting a set of arguments and flags.
% By Chi-Hang Lam, Jan 2019 (based on codes from posting by AJ Johnson)
%
% Usage: 
%   opt = getopt( struct( '<opt>', <default>|'noarg' , '<opt>', <default>|'noarg', ... ), varargin{:} )
% Or
%   [opt vararg] = getopt( struct( '<opt>', <default>|'noarg' , '<opt>', <default>|'noarg', ... ), varargin{:} )
% where
%   varargin = {'<opt>', [<value>], '<opt>', [<value>], ... }
%
% A pair '<opt>', <default> specifies an argument:
%   If '<opt>',<value> is found in varargin,  opt.<opt> = <value>
%   Otherwise, opt.<opt> = <default>. 
%
% A pair '<opt>', 'noarg' specifies a flag:
%   If '<opt>' is found in varargin, opt.<opt> = 1.
%   Otherwise, opt.<opt> = 0.
%
% Options and flags in varargin can appear in arbitrary order.
%
% vararg returns a list of unmatched entries in varargin, 
% which is useful as arguments for a further level of function call.
%  
% EXAMPLE:
% 
% varargin = {'length',10, 'cube', 'density',13}; % typically assigned via a function call instead
% opt = getopt(struct('length',0, 'width',0, 'height',0, 'density',1, 'cube','noarg' ), varargin{:});
% It returns:
%   opt = 
%     length: 10
%      width: 0
%     height: 0
%    density: 13
%       cube: 1
%
% See myplot.m for a further example of usage.
  
outflg = (nargout == 2); % to output unprocessed arg
prop_names = fieldnames(opt);
TargetField = [];
arglist = [];
for ii=1:length(varargin)
  arg = varargin{ii};
  if isempty(TargetField)
    % get field name
    f = find(strcmp(prop_names, arg));
    if ~ischar(arg) | length(f) == 0
      if outflg
        continue; % neglect unknown field
      else
        prop_names
        error('%s ',['invalid property ''',arg,'''; must be one of:'], prop_names{:}); % quit
      end
    end
    TargetField = arg;
    if (ischar(opt.(TargetField))) & strcmp(opt.(TargetField),'noarg')
      opt.(TargetField) = 1;
      TargetField = '';
    end;
  else
    % get value
    opt.(TargetField) = arg;
    TargetField = '';
  end
  arglist = [arglist ii];
end

if ~isempty(TargetField)
  error('Property names and values must be specified in pairs.');
end
for jj=1:length(prop_names) % replace all noarg by 0 for all off-flags
  if ischar(opt.(prop_names{jj})) & strcmp(opt.(prop_names{jj}),'noarg')
    opt.(prop_names{jj}) = 0;
  end;
end;
if outflg
  varargin(arglist) = []; % delete processed arguments
  unprocessed_arg =  varargin;
end;

