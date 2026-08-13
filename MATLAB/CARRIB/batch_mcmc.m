function [] = batch_mcmc(nprocs)
    %%% command: batch_mcmc(5)

    clc;
    close all;
    format compact; format long;
    rng('default');
    warning('off');
    
    if ~exist('./chains', 'dir')
        mkdir('./chains')
    end
    
    %%% load initial parameters
    parameters = load('paramsample.mat');
    parameters = parameters.sample;
    
    %%% options for parallelism
    timeout = 86400;
    nrounds = size(parameters,1); % number of runs
    
    t0 = tic();
    tcluster = parcluster('local');
    tcluster.NumWorkers = nprocs;
    
    %%% number of nodes for each core to work on
    nnodes = zeros(nprocs,1);
    iproc = 1;
    for irounds = 1:nrounds
    nnodes(iproc) = nnodes(iproc) + 1;
    iproc = iproc + 1;
    if( iproc > nprocs )
      iproc = 1;
    end
    end
    
    nstart = zeros(nprocs,1);
    nend = zeros(nprocs,1);
    nstart(1) = 1;
    nend(1) = nnodes(1);
    for iproc = 2:nprocs
    nend(iproc) = nend(iproc-1) + nnodes(iproc);
    nstart(iproc) = nend(iproc-1) + 1;
    end
      
    %%% output information about each rank if we want to
    for iproc = 1:nprocs
    fprintf('Rank %i will work on %i nodes: (%4i - %4i)\n', iproc, nnodes(iproc), nstart(iproc), nend(iproc));
    end
    
    %%% create a job
    fprintf('Creating job ... ');
    j = createJob(tcluster);
    fprintf('done!\n');
    
    %%% create nprocs tasks
    fprintf('Creating %i tasks ... \n', nprocs);
    for iproc = 1:nprocs
        %%% createTask( <job object>, @<function name>, <number of things returned from function>, {list, of, input, arguments})
        fprintf('Creating task %2i/%2i ... ', iproc, nprocs);
        task{iproc} = createTask(j, @mcmc, 1, {nstart(iproc),nend(iproc),parameters});
        fprintf('done!\n');
    end
    fprintf('Done creating tasks\n');
    
    %%% actually submit the jobs
    submit(j);
    
    %%% we can (optionally) wait on these jobs
    % wait(j,'finished');
    wait(j,'finished',timeout); 
    
    %%% now that we know the job is done, fetch the outputs from each task
    results = fetchOutputs(j);
    
    tf = toc(t0);
    fprintf('Integration nodes ...... %.2e\n', nrounds);
    fprintf('Time taken ............. %f seconds\n', tf);
    
end