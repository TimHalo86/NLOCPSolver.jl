
function FindClosestPoint(CurPose, Centerline)
    EuclideanPoint = Centerline .- CurPose
    EuclideanDistance = EuclideanPoint[:, 1].^2 + EuclideanPoint[:, 2].^2
    NearistIdx = argmin(EuclideanDistance)
    distance = sqrt(EuclideanDistance[NearistIdx])
    return NearistIdx, distance
end

function FitPoly33(x, y, z)
    if (size(x, 1)) < size(x, 2)
        x = transpose(x)
    end
    if (size(y, 1)) < size(y, 2)
        y = transpose(y)
    end
    if (size(z, 1)) < size(z, 2)
        z = transpose(z)
    end
    A = [ones(size(x, 1))  x  y  x.^2  x .* y  y.^2  x.^3  x.^2 .* y  x .* (y.^2)  y.^3]
    p = pinv(A) * z
    return p
end

function FindCenterLine(init_pos, center_line, lane_yaw) #init_pos = [x_s; y_s]'
    cline_begin_idx, ~ = FindClosestPoint(init_pos, center_line)
    cline_end_idx = cline_begin_idx + 250#Int(ceil( 1.2*40.0 * 4.0 / 1.0 )   )
    mpc_track_info = zeros(cline_end_idx - cline_begin_idx + 1, 4 ) #x,y,yaw,length
    loop_num = size(center_line, 1)
    for mpc_track_info_idx = 1:1:size(mpc_track_info, 1)
        track_idx = mpc_track_info_idx + cline_begin_idx
        mpc_track_info[mpc_track_info_idx, 1] = center_line[mod((track_idx-1), loop_num)+1, 1]
        mpc_track_info[mpc_track_info_idx, 2] = center_line[mod((track_idx-1), loop_num)+1, 2]
        mpc_track_info[mpc_track_info_idx, 3] = lane_yaw[mod((track_idx-1), loop_num)+1, 1]
        if mpc_track_info_idx == 1
            mpc_track_info[mpc_track_info_idx, 4] = 0
        else
            mpc_track_info[mpc_track_info_idx, 4] = mpc_track_info[mpc_track_info_idx-1, 4] + norm(center_line[track_idx%loop_num+1, 1:2] - center_line[(track_idx-1)%loop_num + 1, 1:2] )
        end
    end
    return mpc_track_info
end

function GetLengthParams(mpc_track_info)
    track_width_candidate = collect(range(-3.5, 3.5, length=20))
    fitting_num = size(track_width_candidate, 1) * size(mpc_track_info, 1)
    lane_length_fitting_x = zeros(fitting_num)
    lane_length_fitting_y = zeros(fitting_num)
    lane_length_fitting_s = zeros(fitting_num)

    data_count = 1
    for length_idx = 1:1:size(mpc_track_info, 1)
        for width_idx = 1:1:size(track_width_candidate, 1)
            lane_length_fitting_x[Int(data_count)] = mpc_track_info[length_idx, 1] + track_width_candidate[width_idx] * cos(mpc_track_info[length_idx, 3] + pi/2)
            lane_length_fitting_y[Int(data_count)] = mpc_track_info[length_idx, 2] + track_width_candidate[width_idx] * sin(mpc_track_info[length_idx, 3] + pi/2)
            lane_length_fitting_s[Int(data_count)] = mpc_track_info[end, 4] - mpc_track_info[length_idx, 4]
            data_count += 1
        end
    end
    cost_to_go_para = FitPoly33(lane_length_fitting_x,lane_length_fitting_y,lane_length_fitting_s)
    return cost_to_go_para
end