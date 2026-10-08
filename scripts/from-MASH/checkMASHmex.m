function checkMASHmex()

% compile PH training algorithm
cd(fileparts(which('trainPH.c')));
try
    mex -R2018a -O GCC='/usr/bin/gcc' trainPH.c vectop.c
catch err
    dispMatlabErr(err);
end

% compile Baum-Welch algorithm
cd(fileparts(which('baumwelch.c')));
try
    mex -R2018a -O GCC='/usr/bin/gcc' baumwelch.c vectop.c fwdbwd.c
catch err
    dispMatlabErr(err);
end

% compile model error calculation
cd(fileparts(which('calcmdlconfiv.c')));
try
    mex -R2018a -O GCC='/usr/bin/gcc' calcmdlconfiv.c vectop.c fwdbwd.c
catch err
    dispMatlabErr(err);
end
    