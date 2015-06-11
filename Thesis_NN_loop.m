% 
% Aron Martinez
%
% Pattern Recognition using Neural Network 
%   Allows setting of desired parameters for NN training
%

%% Intitialization 
clc
clear

% ThesisStart % Uncomment when running from command line

%% Define the inputs and targets
Thesis_main
inputs = nninputs;
targets = nntargets;
size(inputs)
size(targets)

hiddenNeuronsStart = 2; % number of hidden neurons min
hiddenNeuronsEnd = 20;   % number of hidden neurons max
desiredPerc = 85;  % sets the desired train percentage to get correct
trainSets = 30;  % sets number of training session per each hidden layer #

%% Perform loops
for k = hiddenNeuronsStart:hiddenNeuronsEnd
    for i = 1:trainSets
        % Create a Pattern Recognition Network
        hiddenLayerSize = k;
        net = patternnet(hiddenLayerSize);

        % Choose Input and Output Pre/Post-Processing Functions
        % For a list of all processing functions type: help nnprocess
        net.inputs{1}.processFcns = {'removeconstantrows','mapminmax'};
        net.outputs{2}.processFcns = {'removeconstantrows','mapminmax'};

        % Setup Division of Data for Training, Validation, Testing
        % For a list of all data division functions type: help nndivide
        net.divideFcn = 'dividerand';  % Divide data randomly
        net.divideMode = 'sample';  % Divide up every sample
        
        net.divideParam.trainRatio = 70/100;
        net.divideParam.valRatio = 15/100;
        net.divideParam.testRatio = 15/100;

        % For help on training function 'trainscg' type: help trainscg
        % For a list of all training functions type: help nntrain
        net.trainFcn = 'trainscg';  % Scaled conjugate gradient

        % Choose a Performance Function
        % For a list of all performance functions type: help nnperformance
        net.performFcn = 'mse';  % Mean squared error

        % Choose Plot Functions
        % For a list of all plot functions type: help nnplot
        net.plotFcns = {'plotperform','plottrainstate','ploterrhist', ...
          'plotregression', 'plotfit'};

        net.trainParam.epochs = 1000; 	% Sets max # of epochs to ...

        % Train the Network 
        [net,tr] = train(net,inputs,targets);

        % Test the Network
        outputs = net(inputs);
        errors = gsubtract(targets,outputs);
        performance = perform(net,targets,outputs);

        % Recalculate Training, Validation and Test Performance
        trainTargets = targets .* tr.trainMask{1};
        valTargets = targets  .* tr.valMask{1};
        testTargets = targets  .* tr.testMask{1};
        
        trainPerformance = perform(net,trainTargets,outputs);
        valPerformance = perform(net,valTargets,outputs);
        testPerformance = perform(net,testTargets,outputs);

        % Determine value of % error in training
        t = testTargets;
        y = outputs;
        known = find(~isnan(sum(t,1)));

        t = t(:,known);
        y = y(:,known);
        numSamples = size(t,2);
        [c,cm] = confusion(t,y);

        perc = sum(diag(cm))/numSamples*100;

        %% Saving the values of the weights and biases
        if perc >= desiredPerc
            N = hiddenLayerSize;
            weights_I = net.IW{1,1};
            weights_L = net.LW{2,1};
            bias_1 = net.b{1};
            bias_2 = net.b{2};
            testCorrectPerc = perc
            percentage(i,k - hiddenNeuronsStart + 1) = perc; %trainCorrectPerc
        else
            percentage(i,k  - hiddenNeuronsStart + 1) = perc;
        end
        
        i = i+1;
    end
    % max(percentage)
end

percentage
maxPerNN = max(percentage)
avgPerNN = mean(percentage)

save('temp', 'percentage', 'maxPerNN', 'avgPerNN')


