function r = randi_u64
r_u32 = randi(intmax('uint32'),2,1,'uint32');
r = (uint64(intmax('uint32')) + 1) * uint64(r_u32(1)) + uint64(r_u32(2));