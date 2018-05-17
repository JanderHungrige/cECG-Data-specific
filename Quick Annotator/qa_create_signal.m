function result = qa_create_signal(name, unit, sample_time, values)

  % Make row or column vector
  to_row    = @(x) x(:).';
  to_column = @(x) x(:);

  % Verify signal name
  if (isa(name, 'char') && isvector(name))
    result.name = to_row(name);
  else
    error('Name must be a 1D vector of class char.');
  end;
  
  % Verify signal unit
  if (isa(unit, 'char') && isvector(unit))
    result.unit = to_row(unit);
  else
    error('Unit must be a 1D vector of class char.');
  end;
  
  % Verify start time (FOR FUTURE USE)
  % if (isa(start_time, 'char') && is_valid_start_time(start_time))
  %   result.start_time = to_row(start_time);
  % else
  %   error('Invalid start time.');
  % end;
  
  % Verify sample time
  if (isa(sample_time, 'double') && isscalar(sample_time))
    result.sample_time = sample_time;
  else
    error('Sample time must be a 1D vector of class char.');
  end;
  
  % Verify signal values
  if (isa(values, 'double') && isvector(values))
    result.values = to_column(values);
  else
    error('Values must be a 1D vector of class double.');
  end;
  
end   % end function

% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%   Validate time string
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
function tf = is_valid_start_time(utc_time)

  % Examples
  %   2017-03-07T12:40:34+00:00
  %   2017-03-07T12:40:34Z

  % Expressions
  expr_date   = '\d{4}-\d{2}-\d{2}';
  expr_time   = '\d{2}:\d{2}:\d{2}(\.\d{1,8})?';
  expr_offset = '(Z|([+-]\d{2}:\d{2}))';
  
  expr = [ '^' expr_date 'T' expr_time expr_offset '$' ];

  % Test the string
  tf = (regexp(utc_time, expr, 'once') == 1);

end   % end function
