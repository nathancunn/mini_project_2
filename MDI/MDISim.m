function MDISim(uniqueID)


uniqueIdentifier = 2;
nComponents = 50;
nSamples    = 100;
fileName1   = 'Data/MDItestdata1.csv';
fileName2   = 'Data/MDItestdata2.csv';
fileName3   = 'Data/MDItestdata3.csv';
fileName4   = 'Data/MDItestdata4.csv';
fileName5   = 'Data/MDItestdata5.csv';
fileName6   = 'Data/MDItestdata6.csv';
dataType1   = 'TimeCourseUnequalSpacing';
dataType2   = 'TimeCourseUnequalSpacing';
dataType3   = 'TimeCourseUnequalSpacing';
dataType4   = 'TimeCourseUnequalSpacing';
dataType5   = 'TimeCourseUnequalSpacing';
dataType6   = 'TimeCourseUnequalSpacing';
hyperSamplingFrequency = 1;
samplingFrequency = 10;
initialise = true;


MDI('', nComponents, nSamples, {fileName1, fileName2, fileName3,...
    fileName4, fileName5, fileName6},...
    {dataType1, dataType2, dataType3, dataType4, dataType5, dataType6},...
    hyperSamplingFrequency,...
    samplingFrequency, uniqueIdentifier, initialise);




end
