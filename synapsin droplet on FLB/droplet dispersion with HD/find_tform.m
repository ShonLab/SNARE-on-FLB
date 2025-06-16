function [tform, img_comp] = find_tform(img,disp_range)

% Split the input image into two halves
I1 = img(1:end/2, :);
I2 = img(1+end/2:end, :);

% Initialize transformation parameters
Scale = 1;
RotationAngle = 0;
Translation = [0, 0];

% Create a figure and set the callback function
f = figure('KeyPressFcn', @figureCallback);

% Update the display with initial parameters
updateDisplay(I1, I2, Scale, RotationAngle, Translation);

    % Define the figure callback function
    function figureCallback(~, key)
    
        switch key.Key
            case 'a'
                Scale = Scale + 0.001;
            case 's'
                Scale = Scale - 0.001;
            case 'leftarrow'
                Translation(1) = Translation(1) - 0.5;
            case 'rightarrow'
                Translation(1) = Translation(1) + 0.5;
            case 'downarrow'
                Translation(2) = Translation(2) + 0.5;
            case 'uparrow'
                Translation(2) = Translation(2) - 0.5;
            case 'd'
                RotationAngle = RotationAngle - 0.1;
            case 'f'
                RotationAngle = RotationAngle + 0.1;
            case 'escape'
                close(f);
                return
        end
    
        % Update the display with the new transformation parameters
        updateDisplay(I1, I2, Scale, RotationAngle, Translation);
    
    end

    % Define the updateDisplay function
    function updateDisplay(I1, I2, Scale, RotationAngle, Translation)
    
        % Apply the transformation to I2
        tform = simtform2d(Scale, RotationAngle, Translation);
        I1_moved = imwarp(I1, tform, 'OutputView', imref2d(size(I2)));
        
        % Display the images
        img_comp = double(cat(3, I2 / disp_range(2), I1_moved / disp_range(1), zeros(size(I2))));
        imshow(img_comp);  
        title(sprintf('Translation: [%.2f, %.2f], Rotation: %.2f, Scale: %.2f', ...
            tform.Translation(1), tform.Translation(2), tform.RotationAngle, tform.Scale));
    end

   % Wait until the user exits the interactive alignment process
    while true
        pause(0.1); % Check for key presses every 0.1 seconds
        if ~ishandle(f) % Check if the figure window has been closed
            break;
        end
    end

    % Return the transformation matrix and composite image
    tform = simtform2d(Scale, RotationAngle, Translation);
    I1 = imwarp(I1, tform, 'OutputView', imref2d(size(I2)));
    img_comp = double(cat(3, I2, I1, zeros(size(I2))));
end
