classdef Windows < handle
    properties
        idx
        w_size
        ratio
        pre_sig
        post_sig
        final
    end
    methods
        function obj = Windows(pre, post, ws, ratio)
            obj.idx = 1;
            obj.pre_sig = pre;
            obj.post_sig = post;
            obj.w_size = ws;
            obj.ratio = ratio;
            obj.final = false;
            if ratio >= 1
                error('overlap ratio > 1');
            end
        end
        function [pre_frame, post_frame, final] = Next(obj) % Get next frame
            len = length(obj.pre_sig);
            jump = round((1 - obj.ratio) * obj.w_size);
            if obj.final
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
            obj.final = final;
            obj.idx = obj.idx + jump;
        end
    end
end