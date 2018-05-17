% QA_FROM_UNIX_TIME_SECONDS   Convert unix epoch time to a ISO 8601 formatted string.
%
%   Converts a given number of seconds elapsed since 1-jan-1970 (unix epoch)
%   to a ISO 8601 formatted string.
%
function str = qa_from_unix_time_seconds(unix_epoch)

  % Validate argument
  if ~(isa(unix_epoch, 'int64') && isscalar(unix_epoch))
    error('Argument must be a scalar of class int64.');
  end;

  % Invoke external utility
  [ path, ~, ~ ] = fileparts(mfilename('fullpath'));  
  [status, cmdout] = dos(sprintf('%s\\TimeConverter.exe %d', path, unix_epoch));
  
  % Check results
  if (status == 0)
    str = cmdout;
  else
    error('Failed to convert time representation.');
  end;
  
  % -------------------------------------8<-------------------------------------
  
  % Direct invokation of .NET assemblies do not work (yet) as MATLAB fails
  % to load .NET runtime v4.6.2.
  
  % Convert
	%str = System.DateTimeOffset ...
  %        .FromUnixTimeSeconds(unix_epoch) ...
  %        .ToString('O');
  
end   % end function
