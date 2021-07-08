% choose units for labels for rdtx
function units=rdtx_units(name)

switch name
    case 'x'
        units = '\omega_0/c';
    case 'p'
        units = '/ mc';
end