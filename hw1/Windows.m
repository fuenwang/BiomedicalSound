classdef Windows < handle
    properties
        idx
        w_size
        pre_sig
        post_sig
    end
    methods
        function obj = Windows(pre, post, ws)
            obj.idx = 1;
            obj.pre_sig = pre;
            obj.post_sig = post;
            obj.w_size = ws;
        end
        function [pre_frame, post_frame, final] = Next(obj) % Get next frame
            len = length(obj.pre_sig);
            if obj.idx > len
                error('Frame out of bound');
            end
            if obj.idx + obj.w_size - 1 > len
                pre_frame = obj.pre_sig(obj.idx:end);
                post_frame = obj.post_sig(obj.idx:end);
                final = true;
            else
                pre_frame = obj.pre_sig(obj.idx:(obj.idx + obj.w_size - 1));
                post_frame = obj.post_sig(obj.idx:(obj.idx + obj.w_size - 1));
                final = false;
            end
            obj.idx = obj.idx + obj.w_size;
        end
    end
end