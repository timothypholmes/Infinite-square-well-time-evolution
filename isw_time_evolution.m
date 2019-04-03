%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%                 Timothy Holmes
%        infinite square well time evolution
%                   (4/27/18)       
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%isw_time_evolution
%Time evolution for a continuous system

%% Inputs
% m = mass of particle (an electron in this case)
% num = a scalar that specifies the largest energy eigenstate
% x = a vector that sets the range of the well
% psi0 = a column vector that represents the value of the wave function
% at psi(x,0)
%% Outputs
% A video file as a .avi
% An updated plot through time
% Number of iterations in the for loop
%
% For more detailed description of the code
% please see attached pdf

function isw_time_evolution(m,num,x,psi0)

load('data.mat');

%FileData = load('data.mat');
%csvwrite('data.csv', FileData.M);


%% Constants

hbar = 6.58211951*10^-16;
L = x(2001);
n = 1:num;

%% Normalize

A = 1/sqrt(trapz(x,conj(psi0).*psi0)); %normalizing wave function with integration
psi0Norm = A*psi0; %Normalizing wave function

%% Plot
fig = figure;
hold on
plotReal = plot(x,real(psi0), 'linewidth', 2); %real number psi part of the plot
plotImag = plot(x,imag(psi0), 'linewidth', 2); %imag number psi part of the plot
plotAbs = plot(x,abs(psi0), 'linewidth', 2); %abs value of psi part of the plot
legend('real','imaginary','Absolute Value') %labels
xlabel('x (nm)')
ylabel('\bf{\psi(t)}')
ylim([-A A]);

video = VideoWriter('TimeEvolution.avi'); %starts writting frames
open(video);

count = 0;
dt = 2000; %time step
timeTotal = 1*10^8; %total time

%% preallocate 

phi = zeros(length(x),num);
c = zeros(1,length(n));
En = zeros(1,length(n));

%% Time Evolution, eigenstates, c terms, and updates plot

for j = 1:dt:timeTotal
    
    time = j;%*1^-18;
    psi = zeros(size(psi0));
    
    fprintf('Size of psi')
    size(psi)
    
    for k = 1:length(n)
    
        phi(:,k) = sqrt(2/L)*sin((n(k)*pi.*x)/L);
        En(:,k) = (n(k).^2*pi^2*hbar^2)/(2*m*(L^2));
        c(k) = trapz(x,conj(phi(:,k)).*psi0Norm);
    
        psi = psi + c(k).*phi(:,k)*exp(-1i*En(k)*time/hbar);
    
    end

    count = count + 1; %checking iteration number
      
    title(sprintf('Time Evolution Time Frame Number = %g', count))
    
    set(plotReal, 'YData', real(psi)) %update plot
    set(plotImag, 'YData', imag(psi))
    set(plotAbs, 'YData', abs(psi))
    
    currentFrame = getframe(gcf); %grabs frams from iteration
    writeVideo(video,currentFrame); %writes frames out to .avi file
    
    drawnow %draws updated plot
    pause(0.005) %waits for next iteration
                
end


close(fig);
close(video);

end




