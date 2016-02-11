%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pattern Recognition using Neural Network 
%   Allows setting of desired parameters for NN training
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define inputs and targets
load nninoutenergy.mat;
inputs = nninputs; sizeInputs = size(inputs);
targets = nntargets; sizeTargets = size(targets);

%% Initialize all variables
hiddenNeuronsStart = 2; % number of hidden neurons min
hiddenNeuronsEnd = 30;   % number of hidden neurons max
desiredPerc = 0.80;  % sets the desired train percentage to get correct
trainSets = 100;  % sets number of training session per each hidden layer #

conditionmet = 0; % flag for keeping track if desired percentage met 

%% Print conditions for NN training
fprintf('********************************************************\n');
fprintf('* Starting NN training with %d training sessions\n',trainSets);
fprintf('*    from %d to %d number of hidden NN\n', hiddenNeuronsStart, hiddenNeuronsEnd);
fprintf('* Looking for %.2f percentage correct\n', desiredPerc* 100);
fprintf('* Input size is: [%d, %d], target size is: [%d, %d]\n', sizeInputs, sizeTargets);
fprintf('********************************************************\n');

%% Perform loops
for k = hiddenNeuronsStart:hiddenNeuronsEnd
    %fprintf('\n*****\nNumber of hidden neurons = %d\n*****\n', k);

    for i = 1:trainSets
        % Create a Pattern Recognition Network
        hiddenLayerSize = k;
        net = patternnet(hiddenLayerSize);

        % Choose Input and Output Pre/Post-Processing Functions
        net.inputs{1}.processFcns = {'removeconstantrows','mapminmax'};
        net.outputs{2}.processFcns = {'removeconstantrows','mapminmax'};

        % Setup Division of Data for Training, Validation, Testing
        net.divideFcn = 'dividerand';  % Divide data randomly
        net.divideMode = 'sample';  % Divide up every sample
        
        net.divideParam.trainRatio = 90/100;
        net.divideParam.valRatio = 5/100;
        net.divideParam.testRatio = 5/100;

        % Choose the Training Function
        net.trainFcn = 'trainscg';  % Scaled conjugate gradient

        % Choose a Performance Function
        %net.performFcn = 'mse';  % Mean squared error
        net.performFcn = 'crossentropy';  % Cross-entropy

        % Choose Plot Functions
        net.plotFcns = {'plotperform','plottrainstate','ploterrhist', ...
          'plotregression', 'plotfit'};

        % Set the Training Parameters
        net.trainParam.epochs = 1000;
        net.trainParam.showWindow = false;
        net.trainParam.goal = .01;

        % Train the Network 
        net = init(net);
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

        % fprintf('NN Percent: %.2f\n',testPerformance);
        
        %% Save the values of the weights and biases if met
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
            percentage(i,k - hiddenNeuronsStart + 1) = testPerformance;
        end
        
        i = i+1;
    end
end

%% Show Results of NN training
fprintf('Results of NN training:\n');
percentage  % Output all percentages
maxPerNNsort = sort(percentage,'descend');
top5max = maxPerNNsort(1:5,:) % Output top 5 percentages 

% Find max percentage with index to determine # of hidden neurons
maxPerNN = max(percentage);
[maxPercentage, NeuronIndex] = max(maxPerNN);
MaximumPercentage = maxPercentage * 100;
HiddenNeurons = NeuronIndex + hiddenNeuronsStart - 1;

% Print the max test percentage and the number of hidden neurons 
fprintf('Max test percentage is %.2f at %d Hidden Neurons \n',MaximumPercentage, HiddenNeurons);

%% Save desired values
    % Will only save the weights and biases of the last set that met the condition and not all sets
if conditionmet
    fprintf('Optimum hidden layer size = %d', N)
    save('WandB', 'weights_L', 'weights_I', 'bias_2', 'bias_1', 'testCorrectPerc', 'N')
end

save('percPerNN', 'percentage', 'top5max', 'maxPercentage')




