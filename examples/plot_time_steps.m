% Specify the file path
file_path = 'output.txt'; % Replace 'your_file.txt' with the actual file path

% Load the numeric data from the file
try
    data = load(file_path);
    disp('File loaded successfully.');
    
    % Display the loaded data
    disp('Loaded data:');
    disp(data);
    
    % Now you can perform further operations with the loaded data
    
catch
    disp('Error loading the file. Please check the file path and format.');
end

dt1=0.00404858;
dt2=0.00121457;

N=100;
xx=0:1:N;
data(1:N)
figure(1)
plot(dt1*xx,dt1*xx,'bo'), hold on
plot(dt2*xx,dt2*xx,'ko'), hold on
plot(data(1:N+1),data(1:N+1),'r*'), hold off
axis([0,data(N+1),0,data(N+1)])