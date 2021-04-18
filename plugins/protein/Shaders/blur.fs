
//   <namespace name="postprocessing">
//     <shader name="blur">
//       <snippet type="version">330</snippet>
//       <snippet type="string">
//         <!--
//         uniform sampler2D screenTexture;

//         void main()
//         {
//             col = 4.0 / 16 * texelFetch(screenTexture, ivec2(gl_FragCoord.xy), 0).xyz
//             col += 2.0/16 *  texelFetch(screenTexture, ivec2(gl_FragCoord.xy + vec2(-1,0)), 0).xyz;
//             col += 2.0/16 *  texelFetch(screenTexture, ivec2(gl_FragCoord.xy + vec2(1,0)), 0).xyz;
//             col += 2.0/16 *  texelFetch(screenTexture, ivec2(gl_FragCoord.xy + vec2(0,1)), 0).xyz;
//             col += 2.0/16 *  texelFetch(screenTexture, ivec2(gl_FragCoord.xy + vec2(0,-1)), 0).xyz;
//             col += 1.0/16 *  texelFetch(screenTexture, ivec2(gl_FragCoord.xy + vec2(-1,-1)), 0).xyz;
//             col += 1.0/16 *  texelFetch(screenTexture, ivec2(gl_FragCoord.xy + vec2(1,-1)), 0).xyz;
//             col += 1.0/16 *  texelFetch(screenTexture, ivec2(gl_FragCoord.xy + vec2(-1,1)), 0).xyz;
//             col += 1.0/16 *  texelFetch(screenTexture, ivec2(gl_FragCoord.xy + vec2(1,1)), 0).xyz;
//             FragColor = vec4(col, 1.0);
//         }

//       -->
//     </snippet>