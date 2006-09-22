function process_optional_args(varargin)
% PROCESS_OPTIONAL_ARGS  Assign name, value pairs in the callers workspace
%
% This function ensures that if any arguments are passed, that they are passed
% as string, value pairs, then assigns each value to the variable referred to
% in the preceeding argument.
%

% CVS INFO %
%%%%%%%%%%%%
% $Id: process_optional_args.m,v 1.2 2006-09-22 18:00:43 scottl Exp $
%
% REVISION HISTORY
% $Log: process_optional_args.m,v $
% Revision 1.2  2006-09-22 18:00:43  scottl
% updated MSGID for warning message.
%
% Revision 1.1  2006-09-20 21:58:30  scottl
% initial checkin
%

% LOCAL VARS %
%%%%%%%%%%%%%%


% CODE START %
%%%%%%%%%%%%%%
arglen = length(varargin);
if rem(arglen, 2) ~= 0
    error('optional arguments must be variable name and value pairs');
end
for ii=1:2:arglen-1
    assignin('caller', varargin{ii}, varargin{ii+1});
    warning('MBOCR:override', 'overriding value for variable: %s', ...
            varargin{ii});
end
