#version 430
out vec4 fragPosition;

uniform int level;
uniform int level_max;
uniform float weightFactor = 1.0;

uniform sampler2D inputTex_fragPosition;
uniform sampler2D pyramid_fragPosition;
uniform sampler2D maxDepth_texture;


void pull(out vec4 pulledFragPosition) {
    // read all finer pixels' values
    vec4 fp_bl = (texelFetch(pyramid_fragPosition, ivec2(gl_FragCoord.xy) * 2 + ivec2(0), level - 1));
    vec4 fp_tl = (texelFetch(pyramid_fragPosition, ivec2(gl_FragCoord.xy) * 2 + ivec2(0, 1), level - 1));
    vec4 fp_tr = (texelFetch(pyramid_fragPosition, ivec2(gl_FragCoord.xy) * 2 + ivec2(1), level - 1));
    vec4 fp_br = (texelFetch(pyramid_fragPosition, ivec2(gl_FragCoord.xy) * 2 + ivec2(1, 0), level - 1));

    vec4 distanceWeights = vec4(1);
    float distanceWeightsSum = 4.0;

    pulledFragPosition.xyz = fp_bl.xyz * fp_bl.w * distanceWeights.x + fp_tl.xyz * fp_tl.w * distanceWeights.y +
                             fp_tr.xyz * fp_tr.w * distanceWeights.z + fp_br.xyz * fp_br.w * distanceWeights.w;

    pulledFragPosition.w = fp_bl.w * distanceWeights.x + fp_tl.w * distanceWeights.y + fp_tr.w * distanceWeights.z +
                           fp_br.w * distanceWeights.w;

    if (any(isnan(pulledFragPosition.xyz)) || any(isinf(pulledFragPosition.xyz))) {
        pulledFragPosition = vec4(0);
        return;
    }

    if (pulledFragPosition.xyz != vec3(0)) {
        pulledFragPosition.xyz = pulledFragPosition.xyz / distanceWeightsSum;
        pulledFragPosition.w = min(1, pulledFragPosition.w);
    } else {
        pulledFragPosition.w = 0;
    }
}

void main() {
    fragPosition = vec4(0);

    if (level == 0) {
        fragPosition = texelFetch(inputTex_fragPosition, ivec2(gl_FragCoord.xy), 0);
        float fragZ = abs(fragPosition.z);

        float distanceWeight = fragZ * fragZ;

        float maxDepth = texelFetch(maxDepth_texture, ivec2(0, 0), level_max - 1).x;
        if (fragPosition.xyz != vec3(0)) {
            fragPosition.w =
                weightFactor * fragPosition * distanceWeight / (maxDepth * maxDepth); // weight for this pixel
        } else {
            fragPosition.w = 0; // no weight fofragPositionr this void pixel
        }
        return;
    }

    pull(fragPosition);
}