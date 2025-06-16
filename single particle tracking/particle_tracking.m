function out = particle_tracking(img_filtered, threshold, maxdisp, len)
    % initialize
    out = [];
    pre = [];
    last_id = 0;

    for imageNumber = 1:size(img_filtered, 3)
        % particle detection
        [~, xpos, ypos] = particle_detection(img_filtered(:, :, imageNumber), threshold, 0);
        cur = [xpos ypos repmat(imageNumber, numel(xpos), 1) zeros(numel(xpos), 1)];

        if isempty(pre) && ~isempty(cur)
            % first frame processing
            cur(:, 4) = (last_id + 1):(last_id + numel(xpos));
            last_id = last_id + numel(xpos);
        elseif ~isempty(pre) && ~isempty(cur)
            % distance calculation and ID matching
            for idx = 1:size(cur, 1)
                distances = vecnorm(pre(:, 1:2) - cur(idx, 1:2), 2, 2);
                [min_dist, min_idx] = min(distances);

                if min_dist <= maxdisp
                    cur(idx, 4) = pre(min_idx, 4);
                    pre(min_idx, :) = [];
                else
                    % newly detected particles
                    last_id = last_id + 1;
                    cur(idx, 4) = last_id;
                end
            end
        else
            pre=[];
            continue;
        end

        % assign new ID if there are new particles
%         new_particles = cur(:, 4) == 0;
%         cur(new_particles, 4) = (last_id + 1):(last_id + sum(new_particles));
%         last_id = last_id + sum(new_particles);

        % save result and renew current state
        out = [out; cur];
        pre = cur;
    end

    % deletion of particles that found in less than minimum frames(len)
    if ~isempty(out)
        edge=unique(out(:,4));
        id_counts = histcounts(out(:, 4),[edge-0.5;edge(end)+0.5]);
        valid_ids = find(id_counts >= len);
        out = out(ismember(out(:, 4), valid_ids), :);
    end
end
