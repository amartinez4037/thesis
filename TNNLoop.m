%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pattern Recognition using Neural Network 
%   Allows setting of desired parameters for NN training
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define the inputs and targets
load nninoutenergy.mat;

inputs = nninputs; sizeInputs = size(inputs)
targets = nntargets; sizeTargets = size(targets)
conditionmet = 0;

inputs(1:9,1:5);


%% Initialize all variables
hiddenNeuronsStart = 2; % number of hidden neurons min
hiddenNeuronsEnd = 30;   % number of hidden neurons max
desiredPerc = 0.80;  % sets the desired train percentage to get correct
trainSets = 500;  % sets number of training session per each hidden layer #

fprintf('***************************************************\n');
fprintf('Starting NN training with %d training sessions\n',trainSets);
fprintf('   from %d to %d number of hidden NN\n', hiddenNeuronsStart, hiddenNeuronsEnd);
fprintf('Looking for %.2f percentage correct\n', desiredPerc* 100);
fprintf('***************************************************\n');

%% Perform loops
for k = hiddenNeuronsStart:hiddenNeuronsEnd
    %fprintf('\n*****\nNumber of hidden neurons = %d\n*****\n', k);

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
        %net.performFcn = 'mse';  % Mean squared error
        net.performFcn = 'crossentropy';  % Cross-entropy

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

        tind = vec2ind(targets);
        yind = vec2ind(outputs);
        percErrors = sum(tind ~= yind)/numel(tind);

        % Recalculate Training, Validation and Test Performance
        trainTargets = targets .* tr.trainMask{1};
        valTargets = targets  .* tr.valMask{1};
        testTargets = targets  .* tr.testMask{1};
        
        trainPerformance = perform(net,trainTargets,outputs);
        valPerformance = perform(net,valTargets,outputs);
        testPerformance = perform(net,testTargets,outputs);

        %% Incorrectly based on these % at one point
        % Determine value of % error in training
        %t = testTargets;
        %y = outputs;
        %known = find(~isnan(sum(t,1)));
        %TestOut = outputs(tr.testInd); TestTar = targets(tr.testInd);
        %t = t(:,known); y = y(:,known);
        %numSamples = size(t,2);
        %[c,cm] = confusion(t,y);
        %perc = sum(diag(cm))/numSamples * 100;


        % fprintf('NN Percent: %.2f\n',testPerformance);
        
        %% Saving the values of the weights and biases
        if testPerformance >= desiredPerc
            fprintf('Found a top percentage: %d', testPerformance);
            N = hiddenLayerSize;
            weights_I = net.IW{1,1};
            weights_L = net.LW{2,1};
            bias_1 = net.b{1};
            bias_2 = net.b{2};
            testCorrectPerc = testPerformance
            percentage(i,k - hiddenNeuronsStart + 1) = testPerformance;
            conditionmet = 1;
        else
            percentage(i,k - hiddenNeuronsStart + 1) = testPerformance; % perc w/ other
        end
        
        i = i+1;
    end
end

fprintf('Results of NN training:\n');
percentage  % Output all percentages
maxPerNNsort = sort(percentage,'descend');
top5max = maxPerNNsort(1:5,:) % Output top 5 percentages 
maxPerNN = max(percentage);
[maxPercentage, NeuronIndex] = max(maxPerNN);

MaximumPercentage = maxPercentage * 100;
HiddenNeurons = NeuronIndex + hiddenNeuronsStart - 1;

fprintf('Max test percentage is %.2f at %d Hidden Neurons \n',MaximumPercentage, HiddenNeurons);


%avgPerNN = mean(percentage)

if conditionmet
    N
    save('WandB', 'weights_L', 'weights_I', 'bias_2', 'bias_1', 'testCorrectPerc', 'N')
end

save('percPerNN', 'percentage', 'top5max', 'maxPercentage')




