/***************************************************************************************************
 *  @file      tests/geoxy/secular.cpp
 *  @brief     Tests of the place of an ellipse nearest a point, found by its secular equation:
 *             against mpmath's least distances, and no farther than any dense sample.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

#include "secular.hpp"

#include "geoxy/angle.hpp"
#include "geoxy/at.hpp"
#include "geoxy/ellipse.hpp"
#include "geoxy/matrix.hpp"
#include "geoxy/nearest.hpp"
#include "geoxy/place.hpp"
#include "geoxy/point.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace Geoxy
{
  namespace
  {
    /**
     *  An arc of an ellipse, a point, and the least distance between them: the map's x axis, its y
     *  axis and its origin, the start and the sweep in radians, the point, and the distance, by
     *  mpmath at 60 digits from every root of the feet's quartic in e^(it) and the arc's ends,
     *  rounded once.
     */
    using Known = std::array<double, 11>;

    /**
     *  Ellipses turned and with their parameter turned, eccentric to 1e-3, flat as 1e6 to 1, a hair
     *  from a circle, mirrored, a million off the origin, and sheared; whole and in arcs; the point
     *  outside, inside, a hair from the major axis, a hair from the curve, within 1e-6 of the
     *  center, by the evolute's cusp on the major axis, and a million away.
     */
    constexpr std::array Skewed{
        Known{-0x1.f46f313796c4cp-1, 0x1.cb4f45514fbb0p-6, -0x1.fb25b9c9eaee8p-4,
              -0x1.c2294138a5016p-1, 0x1.607ab87690b80p-7, -0x1.7bf133405c19cp-1,
              -0x1.0df0d648f9d34p+1, 0x1.921fb54442d18p+2, 0x1.a2a1255f6b651p+1,
              -0x1.5fdf7aeaf981ap+0, 0x1.2da9b20712229p+1}, // turned, outside, whole
        Known{0x1.a08b19abb28f2p-2, 0x1.c2e7030b1103dp-2, 0x1.9d328802b1d3ap-2,
              0x1.661e0a0f8d071p-1, 0x1.5479474e8fe8cp-1, 0x1.3a7336c538274p-2,
              -0x1.d4bd8bc798d70p+0, -0x1.911787289fccdp+1, 0x1.9384bebd5c292p-2,
              -0x1.bc2da0c93df24p-6, 0x1.fe82a4dd0abebp-5}, // turned, inside, arc
        Known{-0x1.e923fd687f50ep-1, 0x1.0e2ad4fc41054p-4, 0x1.2360b3a413751p-3,
              -0x1.80d4c4bbdb11dp-1, 0x1.8a2678af069eep-1, 0x1.61378e70d7edep-1,
              0x1.1aadb9fb728c6p+1, 0x1.921fb54442d18p+2, 0x1.8bf4b95f10dbep+0,
              0x1.88892f78c6da2p-2, 0x1.5192ef9d02f8cp-3}, // turned, near axis, whole
        Known{-0x1.468e973b98b30p-1, -0x1.87f956c8df649p-1, 0x1.6deb5239a942ep-4,
              -0x1.7520b4d45edf2p-3, 0x1.e2159af405238p-1, -0x1.24604c2d3ef88p-2,
              -0x1.484f82bb19618p-1, 0x1.2af00f5c76398p+2, 0x1.97906145d345dp-1,
              -0x1.75f25b7348ba2p-1, 0x1.aea1e12a9f460p-30}, // turned, near curve, arc
        Known{-0x1.3b55065015bc3p-3, 0x1.6f6772325cfbcp-1, -0x1.e0284872da2f6p-1,
              -0x1.4a58bd6af36b6p-2, -0x1.8c6fe4d178186p-1, 0x1.b152816dd8d20p-5,
              0x1.b2d83268524dcp+1, 0x1.921fb54442d18p+2, -0x1.8c6fc4060c6bcp-1,
              0x1.b152c2622889cp-5, 0x1.71fc66890e526p-1}, // turned, center, whole
        Known{0x1.2d1062ec5641cp-1, 0x1.c7e01c9e4652bp-2, -0x1.f6c80873999f5p-2,
              0x1.b849f3a231cd3p-1, -0x1.6745504fe3104p-1, 0x1.e0394943ca95cp-2,
              -0x1.c1b6d9ecc1604p+1, -0x1.dbafe3915f290p-1, -0x1.0f24cde77b8e6p-1,
              0x1.a5da9a1cbb12fp-6, 0x1.ebf20116b8d6fp-1}, // turned, cusp, arc
        Known{0x1.494729c10fd77p-1, -0x1.f986409fbfc93p-2, -0x1.8cc537def6d39p-3,
              0x1.44f5662494a46p-1, 0x1.65988a93471fep-1, -0x1.56d658628ba94p-2,
              0x1.0d1358ae5cf96p+1, 0x1.921fb54442d18p+2, -0x1.cd6b092c2279ep+17,
              -0x1.da75ada262a4bp+19, 0x1.e847eaeedd6a5p+19}, // turned, far, whole
        Known{-0x1.c7bb944287a0cp-4, 0x1.f6dfa79fdb4f8p-6, -0x1.edd8d001c3646p-1,
              0x1.e61f934482e46p-3, 0x1.2df962f681c40p-5, -0x1.f73ef0cbefab8p-2,
              0x1.91e6d3594d168p-1, 0x1.921fb54442d18p+2, 0x1.3dcafec270be5p-1,
              -0x1.4553597d2e4dap-1, 0x1.48a8ebb4213fbp-9}, // eccentric, inside, whole
        Known{-0x1.42c733cbf7de0p-6, 0x1.a6f6829e80692p-2, -0x1.7ae78cd1ad3ecp-5,
              0x1.d19276b64ffcdp-1, -0x1.073c5dd3e9808p-2, 0x1.10339a73ef2cap-1,
              0x1.7f29e51684d02p+1, 0x1.921fb54442d18p+2, -0x1.3aa66b8ac592fp-2,
              0x1.8711cab7e647cp+0, 0x1.5d3713a466c56p-26}, // eccentric, near curve, whole
        Known{-0x1.4c16c01d0cf3fp-1, 0x1.2117e4dfe66ffp-2, 0x1.4a6b778bb5f5dp-1,
              -0x1.27452dbf9cd86p-2, -0x1.0eb3ef47a6fe2p-1, 0x1.788dbe12c0b84p-1,
              0x1.018ab9d248280p-4, 0x1.0fe9082d680b6p+2, -0x1.0eb3f6effd07cp-1,
              0x1.788dafd3d2969p-1, 0x1.3cb16ae5b4339p-8}, // eccentric, center, arc
        Known{0x1.1a6d6de478202p-6, -0x1.1a355c74f85eap-4, -0x1.8048e5b5eb662p-3,
              0x1.f597f6406de4fp-1, 0x1.0b54206c669aap-1, -0x1.d341888cd94f0p-2,
              0x1.863a7187304b4p+1, 0x1.921fb54442d18p+2, 0x1.55bfeb1848f01p-2,
              0x1.0d326dd33a2bfp-1, 0x1.ed1b19dc4318ep-17}, // eccentric, cusp, whole
        Known{0x1.37cb1b4b66fb8p-1, -0x1.6f7d093840202p-1, 0x1.276aaabb1fe04p-2,
              -0x1.918e4da583d6dp-3, 0x1.b5582092e54bcp-1, -0x1.c512b32bd847ep-1,
              0x1.73033cee045d8p+1, -0x1.67e1e37871fbdp+1, 0x1.d8e6778b12878p+16,
              0x1.e4b0838519964p+19, 0x1.e847ed640b503p+19}, // eccentric, far, arc
        Known{0x1.1bcc3d70daceap+18, -0x1.6b34cd96e6da5p+19, -0x1.abc3d962f0355p+17,
              0x1.11badab09fbaap+19, -0x1.9baad1037e580p-6, -0x1.836ddec6e0d10p-4,
              0x1.cbddba4f981b0p-2, 0x1.921fb54442d18p+2, 0x1.73079056122acp+20,
              -0x1.d390fb666629cp+17, 0x1.4452ad4cb0812p+20}, // flat 1e6 to 1, outside, whole
        Known{0x1.043a359d424bep+18, 0x1.4524106932b5ep+17, -0x1.892456b1bc76fp+19,
              -0x1.eb350b2abbc9fp+18, -0x1.a2ab6d0c13120p-2, -0x1.029c914bc31d4p-2,
              0x1.d11ac29c92cc0p-2, 0x1.eb736c02ccd70p+0, -0x1.896ef971620e4p+13,
              -0x1.eb918bb5d7732p+12, 0x1.d22d1a9e5fd6fp+16}, // flat 1e6 to 1, inside, arc
        Known{0x1.d3e0a0602cdd4p+17, 0x1.d87ef3a6fe8acp+19, -0x1.28477913ebb47p+14,
              -0x1.2b2ffed289c70p+16, -0x1.adeec8d61fdb4p-2, 0x1.0e1bdc166c50cp-1,
              0x1.ca2eb81cf178cp+0, 0x1.921fb54442d18p+2, 0x1.2fb443feb5139p+17,
              0x1.32b3fa071f8c5p+19, 0x1.865b460616bc6p-1}, // flat 1e6 to 1, near axis, whole
        Known{0x1.dd9319d04c70ep+17, -0x1.d030fb8328778p+19, -0x1.73889c045fda3p+15,
              0x1.69214603da8f9p+17, -0x1.ebc2b0b60c3a6p-1, -0x1.fd90e0cb9e260p-2,
              -0x1.111c68f10b728p-1, 0x1.d90924c016338p+1, -0x1.315a7976aea03p+16,
              0x1.28d0dc0f0c839p+18, 0x1.40f5bbd7cbcc6p+2}, // flat 1e6 to 1, near curve, arc
        Known{0x1.cb458ee967e12p+18, 0x1.beb923b347f35p+16, -0x1.9f2e53d9f46f9p+19,
              -0x1.93d53e0e9d9cap+17, -0x1.c851b7d04045ep-1, 0x1.95f0f6efe23b2p-1,
              -0x1.f834589817fd4p+0, 0x1.921fb54442d18p+2, -0x1.092ca6aad6e3cp+0,
              0x1.2cdde77277069p-2, 0x1.193805a666ec6p-1}, // flat 1e6 to 1, center, whole
        Known{0x1.2d1b982c2204ap+19, 0x1.3636eba31c454p+19, -0x1.3c2e5f5c6d20dp+18,
              -0x1.45bdff9e4e6d7p+18, -0x1.b882024bd39f0p-3, -0x1.03e1b11448d90p-2,
              0x1.68dd7f1a72de0p-3, 0x1.237378e0751d0p+1, -0x1.5415fff21f7aap+19,
              -0x1.5e5f061682a42p+19, 0x1.47074375a8493p+14}, // flat 1e6 to 1, cusp, arc
        Known{-0x1.8709ce6a8d26ap+17, 0x1.d467e6bca4690p+19, -0x1.3de812e460deep+15,
              0x1.7ccb8be8360eap+17, 0x1.26172456a58e8p-3, 0x1.90bc803ac13e8p-1,
              -0x1.5aa709cb79848p-1, 0x1.921fb54442d18p+2, 0x1.ebac68a2eb0afp+37,
              0x1.c125330626cafp+39, 0x1.d1a92ef4dbd04p+39}, // flat 1e6 to 1, far, whole
        Known{0x1.f926daba023a0p-1, 0x1.4ce7863c151fcp-3, -0x1.4cc6fcdfb9adap-3,
              0x1.f92fc6408df6dp-1, 0x1.024d7990cf8c8p-2, 0x1.f354dd28d3fdcp-2,
              0x1.d2fe1c17564c8p-1, 0x1.999d343852828p-1, -0x1.f4655d803994ep-1,
              0x1.e3192cea1af15p-2, 0x1.57eec6813c9b7p+0}, // near circle, outside, arc
        Known{-0x1.f49b96b40a3f2p-1, 0x1.ad9eda90d1935p-3, -0x1.ad9ea79aff2d7p-3,
              -0x1.f49b519b8a8dep-1, -0x1.48c774cd0977cp-2, -0x1.1d2619f93e078p-2,
              -0x1.801399e486068p+1, 0x1.921fb54442d18p+2, -0x1.fa2e6fe8f7272p-1,
              -0x1.2b8f272143e9ep-1, 0x1.0fc2a0e97df33p-2}, // near circle, inside, whole
        Known{0x1.d862c88d690a0p-6, -0x1.ffc97681ac53ep-1, 0x1.ffc977ee89696p-1,
              0x1.d86579a0328c0p-6, -0x1.e7f458d68ca74p-2, 0x1.1b0db3f855f94p-1,
              -0x1.ffc8244b8ee94p+1, 0x1.2ea7f0c775c18p+2, -0x1.479efe4678e9bp-1,
              0x1.94f0376fa6c30p-2, 0x1.8bda5d88949dfp-1}, // near circle, near axis, arc
        Known{0x1.d48f35c821a80p-8, -0x1.fffc0d46804b6p-1, 0x1.fff9bf34bc415p-1,
              0x1.d5dd548709600p-8, -0x1.34fd95442be98p-2, -0x1.ca8b8413d4ba0p-3,
              0x1.d13182f3ef218p-1, 0x1.921fb54442d18p+2, -0x1.d494a2f01bce2p-1,
              -0x1.037d932fc307ep+0, 0x1.1bc4a67babd6bp-19}, // near circle, near curve, whole
        Known{-0x1.3d8b340dac1bcp-1, -0x1.91a254f2ddd05p-1, 0x1.91a2550be0098p-1,
              -0x1.3d8b3408bbd9ep-1, -0x1.efb6531167766p-1, -0x1.164654f0a831ap-1,
              -0x1.ece35d5ab5728p+0, -0x1.48af94ef12e52p+2, -0x1.efb63bd7c2b70p-1,
              -0x1.16466cdd7194ep-1, 0x1.ffffdea0addfep-1}, // near circle, center, arc
        Known{-0x1.f84f47eedf124p-4, -0x1.fc1abc21c2a14p-1, 0x1.fc1abbe6885c3p-1,
              -0x1.f84f4eb0ca6bcp-4, -0x1.aa5e76e67fed0p-1, -0x1.7a85f835858b8p-3,
              0x1.e2c46edd754b0p-1, 0x1.921fb54442d18p+2, -0x1.aa5e75e799090p-1,
              -0x1.7a85fdf96574ap-3, 0x1.fffffe3e1c199p-1}, // near circle, cusp, whole
        Known{-0x1.f29bb1a78e5bfp-2, -0x1.bf341d55afa14p-1, 0x1.bf342ef7292c4p-1,
              -0x1.f29bdfe7da5ebp-2, -0x1.816f56786c36cp-1, -0x1.00a634d8d5b7ap-1,
              0x1.b1d1ad0ff6270p+1, -0x1.47912eb4750c0p-1, -0x1.e146ea2ea2de6p+19,
              0x1.49a35d13056cap+17, 0x1.e84802eb5c29dp+19}, // near circle, far, arc
        Known{0x1.8eea92963e0e4p-1, -0x1.250fbb57e8456p-1, -0x1.c57b9d1db9cdep-3,
              0x1.0a89f780051cap-3, 0x1.c2ad74ac4d3d4p-1, 0x1.7e76d8a3414c8p-3,
              -0x1.09cb3c9b91982p+1, 0x1.921fb54442d18p+2, -0x1.04434093b1928p+0,
              0x1.b034f7ebd275cp-3, 0x1.392aa6542ce39p+0}, // mirrored, outside, whole
        Known{0x1.ed63654485015p-1, -0x1.3493ffa369016p-3, -0x1.de3a2cac12281p-3,
              -0x1.2dea7ae8fc9dap-4, 0x1.a601aa758553cp-1, -0x1.302c33663f24cp-2,
              -0x1.0454e5bdc6f06p+1, 0x1.d71150c1a22f4p+1, 0x1.2d1810c74fda7p+0,
              -0x1.0acfe233ac95fp-2, 0x1.203fe74a44f09p-6}, // mirrored, inside, arc
        Known{-0x1.da169902c3b0dp-4, 0x1.1abc3ecfc8da0p-1, 0x1.8fba81cc09f5fp-2,
              -0x1.7a63d5bf6ce4dp-1, 0x1.cb679575e4cc0p-2, -0x1.b41d1c5ef0a30p-2,
              0x1.8dc1a05729884p+0, 0x1.921fb54442d18p+2, 0x1.2a1c6ddf4da20p-2,
              -0x1.b26500107ea90p-5, 0x1.e62d3c88c0418p-4}, // mirrored, near axis, whole
        Known{0x1.4f58990032262p-4, -0x1.5a9c6a21579b4p-4, -0x1.edcc5fa417c0fp-2,
              -0x1.c02e7868dbabcp-1, 0x1.91fd774df8bcep-1, -0x1.a5e2839919e0cp-1,
              -0x1.e3e739a1588ccp+0, 0x1.921fb54442d18p+2, 0x1.91fd7a1fdb9cbp-1,
              -0x1.a5e29e85846f2p-1, 0x1.ccb076586285bp-4}, // mirrored, center, whole
        Known{0x1.87547a455908ep-2, 0x1.4319af317c5a0p-1, -0x1.43045cb4606ecp-5,
              -0x1.6ffc465ee8304p-1, 0x1.bc9c8ff195e90p-3, 0x1.bbff52a74ec60p-5,
              -0x1.d212656a9c1a2p+1, -0x1.6ddd5370423d2p+2, 0x1.ffbde294dbbaep-2,
              0x1.e57918e7986cdp-1, 0x1.ff1a5f7bfe965p-5}, // mirrored, cusp, arc
        Known{0x1.674de6887010ap-1, -0x1.32e9d7ba77c6dp-1, 0x1.fb1b7a4230166p-3,
              -0x1.34ad62fa88ed5p-2, -0x1.90f3c89203fe0p-2, 0x1.ce01627a78920p-4,
              0x1.8f128e866ba16p+1, 0x1.921fb54442d18p+2, 0x1.169096807608cp+19,
              -0x1.9106174ccdbfep+19, 0x1.e847e0d78f46ep+19}, // mirrored, far, whole
        Known{0x1.b199cc37d5b6ep-2, 0x1.3e95c1f6a0797p-1, -0x1.e8e308742f03dp-2,
              -0x1.db8c54c3cd8e1p-2, 0x1.779a9342a2f0bp+18, -0x1.ddcdced9c8642p+19,
              -0x1.8fa6d49eaaa80p+1, 0x1.91f62fbfba220p+0, 0x1.779b21493842ep+18,
              -0x1.ddcd90290648ap+19, 0x1.240e3748c74f4p+1}, // far off, outside, arc
        Known{-0x1.3489ed65a2d04p-2, -0x1.e6b1a932a07bap-1, 0x1.a615482e21383p-1,
              -0x1.326400ee08db8p-2, 0x1.fa0544d4a7b40p+18, -0x1.9c384704dd7f8p+18,
              -0x1.1c9802662c180p-2, 0x1.921fb54442d18p+2, 0x1.fa055124ba89ep+18,
              -0x1.9c38496981c95p+18, 0x1.5b0c0231e4f31p-1}, // far off, inside, whole
        Known{0x1.77b1191a6e08fp-1, 0x1.28d2b42fdfaeep-1, -0x1.0479b13c25e08p-1,
              0x1.a0fddf8cb225ep-1, 0x1.43f238ddf9060p+18, -0x1.9793bb49ef2f2p+18,
              0x1.33a3b5d5d7106p+1, -0x1.7c1abdd4bf9c8p+0, 0x1.43f2399ee52cbp+18,
              -0x1.9793ad77731c1p+18, 0x1.914253b42ca37p-1}, // far off, near axis, arc
        Known{0x1.4bba8db6cf94ap-1, 0x1.2ccb24a05e465p-2, 0x1.a1288e097c9bep-2,
              0x1.4abe6f7a5531ep-1, 0x1.a65c53e7565a0p+19, -0x1.a182a7cf18e44p+19,
              -0x1.605d2f1cfca20p-3, 0x1.921fb54442d18p+2, 0x1.a65c4c1674d4ap+19,
              -0x1.a182a2040a832p+19, 0x1.09d7a8e11e2ddp-15}, // far off, near curve, whole
        Known{0x1.2249fd250e38fp-2, -0x1.bb1dcecaaa8cbp-3, 0x1.ead47adc9c7d9p-1,
              0x1.3c46c07b863e4p-5, 0x1.6e4d30f22c09bp+18, -0x1.72677b7b69bcfp+19,
              0x1.76bd76d549da0p-3, 0x1.921fb54442d18p+2, 0x1.6e4cf4047cf0bp+18,
              -0x1.72677ab41429bp+19, 0x1.869c69c526f59p-5}, // far off, cusp, whole
        Known{-0x1.8916bf4543080p-3, -0x1.c9626315f6717p-1, 0x1.9343fe3448c5ep-1,
              -0x1.8c86b6633b4adp-2, 0x1.affe44645f940p+17, -0x1.4d99f56cf99acp+17,
              -0x1.c1869f08cc368p+1, -0x1.1caa97c6a8157p+2, -0x1.69e685244a9ecp+19,
              -0x1.b03bec1de7679p+18, 0x1.e847e708028aap+19}, // far off, far, arc
        Known{0x1.d3fd9c9a458ecp-2, -0x1.5317ae71f128bp-2, -0x1.a174b2014a818p-2,
              0x1.7e698e05abc51p-1, 0x1.85e0d70097140p-6, 0x1.c51354ace4986p-1,
              0x1.9ef72c35a0d40p-2, 0x1.921fb54442d18p+2, 0x1.40a5cc7071985p-1,
              0x1.8bb5305c45934p+0, 0x1.586fb02e20cb1p-1}, // turned, outside, whole
        Known{0x1.2c6dae64731edp-1, -0x1.1e5c66953a8f8p-2, 0x1.4d7d1804e528dp-1,
              0x1.70ee1af47ac80p-1, -0x1.abc30ed597490p-2, -0x1.b1fe46f84aad4p-2,
              0x1.18308d9daa2f0p+1, 0x1.5540aeb368cbcp+1, -0x1.51acb2e005185p-1,
              -0x1.f10cda01ae6c5p-1, 0x1.855265d57a4b4p-2}, // turned, inside, arc
        Known{-0x1.667b1380cedb5p-2, 0x1.01d97a26c8595p-2, -0x1.d944db4f052b2p-1,
              0x1.0295a35c64e4fp-4, -0x1.7d44c54069c1ap-1, 0x1.367e8259bf1e4p-2,
              0x1.a119a9b534aaap+1, 0x1.921fb54442d18p+2, -0x1.41c9d09fc6a38p+0,
              0x1.88e88931e0afcp-2, 0x1.6db5fbb1f5131p-3}, // turned, near axis, whole
        Known{-0x1.4322d5861c9ecp-1, 0x1.917eab0bc076ap-2, -0x1.4a6f62d7f94bfp-1,
              -0x1.7cf90724149f0p-1, -0x1.dbcc473e91aacp-2, -0x1.c7002dd8b67a4p-1,
              -0x1.fc62ceda739d0p-2, -0x1.36ee29e23c514p+1, 0x1.7a343a31bcfbep-4,
              -0x1.97d5bad0b04a2p-4, 0x1.58cf850cbd318p-24}, // turned, near curve, arc
        Known{-0x1.eaa759f9e23ecp-1, -0x1.21ab57fb47596p-3, 0x1.1cde2e3a16facp-2,
              -0x1.0dcb7153cbda5p-2, -0x1.145fad59d27ccp-1, 0x1.376626910e2c0p-6,
              0x1.d431fe5981798p-1, 0x1.921fb54442d18p+2, -0x1.145fb01d19d36p-1,
              0x1.37662606e2e2fp-6, 0x1.2ad691f249332p-2}, // turned, center, whole
        Known{0x1.c7f1108b23208p-1, 0x1.b2835d423c04bp-2, -0x1.055de50dbf65ap-2,
              0x1.56b69ffeca4b8p-3, -0x1.f99469f3a2b8ap-1, -0x1.d9ae8909f7560p-4,
              0x1.f7befae215638p-1, 0x1.86f650c680940p+1, -0x1.d8ef6d722a0a7p+0,
              -0x1.eb1764c67bfb7p-2, 0x1.0f390a7f0bcbcp-4}, // turned, cusp, arc
        Known{-0x1.7cb57c163b6c0p-1, 0x1.ce0a7ee434e86p-3, -0x1.458e787a5688bp-1,
              0x1.af577c4fec1d8p-5, -0x1.8c0af6f2a7c4ep-1, 0x1.4259d6afb4ac0p-4,
              -0x1.4244cf2c62fd8p-1, 0x1.921fb54442d18p+2, -0x1.dfb4972f3ac73p+19,
              -0x1.6c740465f4e1fp+17, 0x1.e847e2769f61cp+19}, // turned, far, whole
        Known{0x1.a560965a26d16p-1, -0x1.b63ab36280f47p-2, 0x1.5a312f3f920fap-2,
              -0x1.466c7c6502e62p-3, -0x1.447d3d74dc598p-3, -0x1.ad8d82d2d5eecp-1,
              0x1.31f33f6219910p+0, 0x1.0f091ed0fe590p+2, -0x1.236abce149acap-2,
              -0x1.b98f27ac8b170p-1, 0x1.0a4b24b542335p-4}, // eccentric, outside, arc
        Known{-0x1.65654ad01752fp-1, -0x1.56dd257734f10p-1, -0x1.0169dfa1b5600p-3,
              -0x1.df54e5d30de16p-3, 0x1.631d43bca7914p-2, -0x1.6e95525708f0ap-1,
              -0x1.20018f986fc60p-2, 0x1.921fb54442d18p+2, -0x1.1237e28a0841dp-7,
              -0x1.0d2a625bd6d96p+0, 0x1.c3414684558efp-5}, // eccentric, inside, whole
        Known{-0x1.124d471fdc517p-1, -0x1.a269fd3dc85b5p-1, -0x1.ccc1e673b0107p-4,
              -0x1.713b3fe087e30p-3, -0x1.2275c168dd01ap-1, 0x1.d2d4287171528p-2,
              0x1.4db97b4a2901cp+1, 0x1.921fb54442d18p+2, -0x1.005a1ec13c10fp-1,
              0x1.1939f97f24718p-1, 0x1.b700de7d6849bp-25}, // eccentric, near curve, whole
        Known{-0x1.4005410ff9995p-2, 0x1.6b6b429c6fb62p-2, 0x1.1355d07899230p-1,
              -0x1.65a54afad69cep-1, -0x1.746f28d8d1ce0p-5, 0x1.dae5412d95d9ep-1,
              0x1.135dff5fc9daep+1, 0x1.3438e80943fdcp+2, 0x1.d0d67bec49cd2p+19,
              -0x1.2af5aadbe2502p+18, 0x1.e847e5737c70ep+19}, // eccentric, far, arc
        Known{-0x1.a5ec610187fecp+19, -0x1.e8673820af6bdp+18, 0x1.7ef9b808247b3p+15,
              0x1.bb4cea414550fp+14, -0x1.967cd4067d12cp-1, 0x1.e9b64a9c0d6f6p-1,
              0x1.c90b5758f5b44p+1, 0x1.921fb54442d18p+2, 0x1.2b9b85f346470p+20,
              0x1.89d8b7b7efd52p+20, 0x1.1d8c05f4dbf9ep+20}, // flat 1e6 to 1, outside, whole
        Known{-0x1.5a859e029f4fbp+18, -0x1.6a648edc2b257p+19, -0x1.defc18ced050ap+17,
              -0x1.f4ecdb8fedf6dp+18, -0x1.f0b19eb3046a0p-5, 0x1.2d6d3c2a65d54p-2,
              -0x1.56e53072aefc8p+1, -0x1.7d822ccd93cfep+1, 0x1.91d4016a5c883p+18,
              0x1.a43bb203bae34p+19, 0x1.d1173a8cb18c5p-2}, // flat 1e6 to 1, inside, arc
        Known{0x1.9a3467181441ep+18, -0x1.553b1a62b2988p+19, 0x1.2343f5dc90a5bp+18,
              -0x1.e4943a0b6ecf2p+18, 0x1.5f01902ef071ep-1, -0x1.32020424d3e88p-2,
              -0x1.3ed85ce07cb10p-1, 0x1.921fb54442d18p+2, -0x1.43d0bef27b946p+14,
              0x1.0d5fa7564baaep+15, 0x1.ff95e4c517d1ap-1}, // flat 1e6 to 1, near axis, whole
        Known{0x1.b9b9aba67aa25p+18, -0x1.8eb1ae3fd4256p+18, 0x1.1f69223adfc5dp+19,
              -0x1.036940207e047p+19, -0x1.dd610516f6380p-7, 0x1.e0d310754cb78p-3,
              0x1.c95f178eedb20p+1, 0x1.c291ccf1b5ad0p+1, 0x1.546a5f23231cdp+19,
              -0x1.33413c0e7a1f0p+19, 0x1.83464d0ebaba5p+1}, // flat 1e6 to 1, near curve, arc
        Known{-0x1.6f239e8d51a7fp+18, -0x1.d7956fe8ae12dp+18, 0x1.da72d40cc91afp+18,
              0x1.30b58707dd6e7p+19, 0x1.d7ae8d8678834p-2, 0x1.b4050a9fc6f9ep-1,
              0x1.8a0c8ce989800p+1, 0x1.921fb54442d18p+2, 0x1.68743b0996637p-1,
              0x1.1743f470fd408p+0, 0x1.e8ee8285a2fa5p-1}, // flat 1e6 to 1, center, whole
        Known{0x1.bf69d9794d2e1p+15, -0x1.9fd2a766c330ap+14, 0x1.b9ece811a5de7p+19,
              -0x1.9ab48b6cdd200p+18, -0x1.73152c18df8f8p-2, -0x1.9201bf150d670p-4,
              -0x1.4e534104760f8p+0, 0x1.5d18e752090acp+0, -0x1.bacf2121a7a49p+19,
              0x1.9b86bc94e9097p+18, 0x1.a0915f49ace30p+15}, // flat 1e6 to 1, cusp, arc
        Known{-0x1.052c409204c7ap+19, -0x1.2e1c38a56e57bp+16, -0x1.9699f08a7dfbbp+19,
              -0x1.d6575298d33c7p+16, -0x1.a64d9c2935eacp-1, 0x1.4fd093bb5c3dcp-2,
              -0x1.edd3e66bbd550p+0, 0x1.921fb54442d18p+2, 0x1.e0e61eeced10dp+38,
              -0x1.8ec75c6381948p+39, 0x1.d1a93e44e40aep+39}, // flat 1e6 to 1, far, whole
        Known{-0x1.e3e94c352b550p-3, 0x1.f1808b48e6f82p-1, -0x1.f1808b3abeb69p-1,
              -0x1.e3e94c1b1f95bp-3, -0x1.390850cde3f9ep-1, 0x1.5ffd30b320b68p-3,
              -0x1.1c705bfbc9b4cp+1, 0x1.7c80ffb2178a0p+1, 0x1.09a27bfe81e9ep+0,
              -0x1.59a17a686eed6p-3, 0x1.5e207bc915c02p-1}, // near circle, outside, arc
        Known{0x1.3ef27c5f19206p-2, 0x1.e687df0a25c0ap-1, -0x1.e687df22a2a58p-1,
              0x1.3ef27c4765fbep-2, -0x1.9a8bcce525418p-1, -0x1.0d3d471c717fap-1,
              0x1.d54d940af9c70p-1, 0x1.921fb54442d18p+2, -0x1.671999a1131ccp-1,
              0x1.83b315c171e80p-5, 0x1.ac1bc2bde3355p-2}, // near circle, inside, whole
        Known{-0x1.fdc45869d46efp-1, 0x1.6f230aee969a5p-4, -0x1.6f303f93613f3p-4,
              -0x1.fdf03cf8764a1p-1, 0x1.eab3af6a44e2ap-1, -0x1.58397459837dcp-2,
              -0x1.342ea87c56b50p-2, 0x1.d035f41b167c0p-1, 0x1.f4f8748689858p-1,
              -0x1.62e28e2689561p-6, 0x1.e4b57636bab58p-1}, // near circle, near axis, arc
        Known{0x1.b08be11280c3cp-1, 0x1.11f2bc4bd7f6cp-1, -0x1.11f2bc5850730p-1,
              0x1.b08be117aad65p-1, -0x1.bd01d9a0f81a8p-3, -0x1.f7cc8ce3e9700p-1,
              0x1.5ea872cb26132p+1, 0x1.921fb54442d18p+2, 0x1.84f62f2679424p-1,
              -0x1.8a927e4b4a6dep-1, 0x1.a36d631c07299p-31}, // near circle, near curve, whole
        Known{-0x1.bbb9b285a6218p-1, 0x1.fee3743cac50fp-2, -0x1.fee375035b844p-2,
              -0x1.bbb9b4095fc5dp-1, 0x1.bd2496d8c24b0p-1, -0x1.329bfbea33f7ep-1,
              0x1.ef1e553721c26p+1, 0x1.921fb54442d18p+2, 0x1.bd249808ced62p-1,
              -0x1.329bf905b0e8fp-1, 0x1.fffffcdf2c491p-1}, // near circle, cusp, whole
        Known{-0x1.491533a2b839dp-3, 0x1.5dce7aa9ade8ap-2, 0x1.6b4fee6a8dee8p-1,
              -0x1.3969af529109cp-1, 0x1.2791667a7f7a0p-2, -0x1.725d1d16cd4e0p-1,
              0x1.5966279ea5b20p-1, 0x1.921fb54442d18p+2, 0x1.bc88a7eb92f44p+0,
              0x1.986a44e190a6ap+0, 0x1.4626e29c75d1ap+1}, // mirrored, outside, whole
        Known{-0x1.f97af03609108p-1, 0x1.1886b2b89c477p-3, -0x1.3df74428cfefep-4,
              0x1.db38c66de6f17p-6, -0x1.a3a8b31a66c34p-2, 0x1.730cc8175b35cp-2,
              -0x1.1e6c9c54d6310p-2, -0x1.6c51306031286p+2, -0x1.5f14926870e1ap-2,
              0x1.697b330f51fadp-2, 0x1.258426185a57dp-6}, // mirrored, inside, arc
        Known{0x1.2afeb3c4ff404p-2, 0x1.c5b8fd46e4ec9p-1, -0x1.6a02c3c147479p-4,
              -0x1.65f0d2c88d118p-2, 0x1.7080ec81281aep-1, -0x1.ce003537a3254p-2,
              -0x1.0b673716dbf2ap+1, 0x1.921fb54442d18p+2, 0x1.f96aeee1c73a8p-2,
              -0x1.28c96dd8b3b8cp+0, 0x1.04190d88a7eeap-6}, // mirrored, near axis, whole
        Known{-0x1.7537c58445c38p-4, -0x1.7917441bb3517p-2, 0x1.abc3a792d5c42p-4,
              0x1.d74647c315478p-1, 0x1.7fc45ddac8848p-3, -0x1.c96d59815cad0p-4,
              0x1.553a311f7b134p+0, -0x1.4dd1cdf8f9220p-2, 0x1.6a4f2523ada66p-3,
              -0x1.04ac5a69819a1p-1, 0x1.f5b23d8647418p-1}, // mirrored, near curve, arc
        Known{-0x1.5c885e4f0d5b7p-1, -0x1.0c18c1f8e6fbep-1, -0x1.05b7a4f1b6919p-1,
              -0x1.4a7a0e3951a7cp-3, 0x1.5b6136a3ab2e4p-2, -0x1.ea7f0868fb918p-3,
              0x1.452432e9e7628p+0, 0x1.921fb54442d18p+2, 0x1.5b61794966877p-2,
              -0x1.ea7f5035396ffp-3, 0x1.4333a1c7edeb3p-3}, // mirrored, center, whole
        Known{-0x1.1dba9e8c6de85p-2, 0x1.3b64077dd81f7p-1, -0x1.07dffe3e8935dp-2,
              0x1.61c7452cbdd2bp-1, -0x1.f90fb3144988cp-1, 0x1.03482bfa94d18p-2,
              -0x1.5fe5b9dc1e962p+1, 0x1.921fb54442d18p+2, -0x1.64fa39dcea428p+19,
              0x1.4d2551a697a03p+19, 0x1.e847e2eec2ef3p+19}, // mirrored, far, whole
        Known{0x1.eb7ea4402ad25p-2, -0x1.c11232701f28cp-1, 0x1.5802b57e5018ep-1,
              0x1.860ce0462a677p-2, 0x1.5baafb9317980p+18, -0x1.89dc4eda063bcp+19,
              -0x1.458ca8e40501ap+1, 0x1.59b3650745e68p-1, 0x1.5baa63e1f557ep+18,
              -0x1.89dc3d1270f78p+19, 0x1.91a8e6033426ep+0}, // far off, outside, arc
        Known{0x1.ac5cf0a4f6a2fp-1, -0x1.944b52a5595b3p-2, -0x1.0738b21babda8p-2,
              0x1.50f3106937d06p-2, 0x1.c310e33a12ec0p+19, -0x1.88d62b197c108p+19,
              0x1.cfe003ccf6f50p-1, 0x1.921fb54442d18p+2, 0x1.c310e947e7afbp+19,
              -0x1.88d632d169dcap+19, 0x1.950a1ed748a7dp-5}, // far off, inside, whole
        Known{0x1.42bee2b27b010p-2, 0x1.d1174ef474ce3p-2, -0x1.771f5d1d4659ap-1,
              -0x1.c43e1f849a1f2p-2, 0x1.afc395b6f269cp+19, -0x1.f61f052eae570p+16,
              0x1.34489222120c4p+0, -0x1.65a28cc646d14p+0, 0x1.afc3926d03e70p+19,
              -0x1.f61f19b1c2408p+16, 0x1.88eb30a398acdp-3}, // far off, near axis, arc
        Known{0x1.9f1a4fe00a63bp-2, -0x1.dc1a0d1fe2948p-6, -0x1.a9e4bf04aaba4p-3,
              0x1.f0d2c54ddf21ep-1, 0x1.3d8986aeaae4bp+19, -0x1.2fa43eb095f11p+19,
              0x1.96153b2db7f40p-4, 0x1.921fb54442d18p+2, 0x1.3d89786588005p+19,
              -0x1.2fa42a8ea3f6fp+19, 0x1.250461aa04059p-17}, // far off, near curve, whole
        Known{0x1.ea833cbe6df01p-1, 0x1.b39535750a23dp-3, -0x1.b18c2087214cdp-3,
              0x1.2a66de12ebadep-4, 0x1.2c3788f29afbep+19, -0x1.bc061f308ed62p+18,
              0x1.983f1017e4ac0p+0, -0x1.7835cc9ab1953p+2, 0x1.2c3788f29a3cep+19,
              -0x1.bc061f3092927p+18, 0x1.d64b8a6ad4075p-4}, // far off, center, arc
        Known{0x1.10c194111dbe4p-1, -0x1.958c7a92e1bf9p-1, 0x1.9fed16ba37191p-1,
              0x1.85073386d1e34p-2, 0x1.e89fc38ff7f87p+17, -0x1.a7f6fe4dad3f2p+18,
              0x1.8027865cdc1f0p+1, 0x1.921fb54442d18p+2, 0x1.e89fe42f1e855p+17,
              -0x1.a7f70667087d5p+18, 0x1.6e5001c95441ep-1}, // far off, cusp, whole
        Known{-0x1.008a84b4dbce9p-1, -0x1.89bce0c657c20p-1, 0x1.657cae6d41522p-1,
              -0x1.696e096ac9a88p-4, 0x1.2199f32be3e8ep+19, -0x1.a66d1535c35cap+19,
              0x1.bd89fa6c71bcap+1, 0x1.5078af7241308p+1, 0x1.357795b061cfbp+20,
              -0x1.ef6dd7065b4d8p+16, 0x1.e847e8fd16da3p+19}, // far off, far, arc
        Known{-0x1.6dbcda8f5a59ep-1, -0x1.1f4842bc0ce9bp+0, -0x1.f0f9462673af1p-2,
              -0x1.6dd24bdec09d1p-1, 0x0.0p+0, 0x0.0p+0, -0x1.baaaae8c00c0ap+1,
              0x1.921fb54442d18p+2, 0x1.4e45d239ffbf0p-12, 0x1.00c7940cb0592p-10,
              0x1.5d6b96e81d994p-6}, // sheared, whole
        Known{0x1.72c42556e88bep+0, 0x1.b44dfc1b0c647p-3, 0x1.58574612cdf9ap-2,
              0x1.d22b55f2e7e5ep-5, 0x0.0p+0, 0x0.0p+0, 0x1.3fdc571658158p+1, -0x1.780b76444e65dp+2,
              -0x1.124a30c6ccae6p+0, 0x1.53c3de9d63bbap+0, 0x1.76cefbfcf1ee3p+0}, // sheared, arc
        Known{0x1.53849533798d6p-3, 0x1.21a4e169d66c9p-2, 0x1.7f2644d5e23cep-2,
              0x1.418f805126be5p-1, 0x0.0p+0, 0x0.0p+0, 0x1.9a95760d8aa80p+0, 0x1.921fb54442d18p+2,
              -0x1.7f0ea84af9282p+0, 0x1.ae3ba02db440cp+0, 0x1.1263d4c9d88a2p+1}, // sheared, whole
        Known{-0x1.a64bc384f4137p-2, -0x1.71b2c37c1bf8cp-1, -0x1.41e660db56a8ap-3,
              -0x1.af11162ab3629p+0, 0x0.0p+0, 0x0.0p+0, 0x1.9fadf63f3e6aep+1, 0x1.675f0149e7e20p+0,
              -0x1.58f90e0249c64p-11, 0x1.29de28d2e3074p-12, 0x1.fc234fafc6b29p-1}, // sheared, arc
        Known{-0x1.249a163b22555p-1, -0x1.cd525cb651599p-4, -0x1.1099f474b7e14p+0,
              0x1.adda928062db2p-6, 0x0.0p+0, 0x0.0p+0, 0x1.77c753b3f5dd4p+1, 0x1.921fb54442d18p+2,
              -0x1.9ab8146e5594ep+0, 0x1.e5e7bb84da518p-2, 0x1.45c1ced54d5c1p-1}, // sheared, whole
        Known{-0x1.6b9fc32fea05ep+0, -0x1.86c4970c8d18fp-1, -0x1.b88cec7a9721dp-1,
              -0x1.2ec55c2e0dfc8p+0, 0x0.0p+0, 0x0.0p+0, 0x1.113289f5ea000p+1, 0x1.07faf4a790820p+0,
              0x1.5ae7aff72ddc0p-3, 0x1.eb2f1026718ecp-1, 0x1.137284460a98ep+0}, // sheared, arc
        Known{0x1.e6449ae22e470p-1, 0x1.0d247087407a1p+0, 0x1.649a426a047cbp-1,
              0x1.8695d49507c82p-1, 0x0.0p+0, 0x0.0p+0, -0x1.2fb98c13d6d4ap+1, 0x1.921fb54442d18p+2,
              0x1.d6e72726031f6p-11, 0x1.e2bcc6ae1c6c0p-16, 0x1.eca8298cb99d9p-9}, // sheared, whole
        Known{0x1.97f09001f9b2ep-2, -0x1.2fc9f86904043p-4, 0x1.ecb295297513ep-1,
              0x1.b756c4109614ep-2, 0x0.0p+0, 0x0.0p+0, 0x1.37e1575f8bd9cp+1, -0x1.d45e62d4f389bp+1,
              0x1.958d3d8266d28p+0, -0x1.b5418aaea1b78p-2, 0x1.d70f9f2f37700p-1}, // sheared, arc
        Known{-0x1.eafe988fa1231p-1, -0x1.12dec4041b7bbp-1, -0x1.7326280598952p-2,
              -0x1.94c5d591457f2p-3, 0x0.0p+0, 0x0.0p+0, -0x1.13070d5248690p+0,
              0x1.921fb54442d18p+2, -0x1.8170826f73646p+0, -0x1.29dbee1184db2p+0,
              0x1.8623779859e57p-1}, // sheared, whole
        Known{0x1.d47e1db5589c8p-1, -0x1.c5afeca127604p-2, 0x1.f06e046e36f0dp-2,
              -0x1.003d527660f82p-2, 0x0.0p+0, 0x0.0p+0, 0x1.bd1ce72b99c40p-3, 0x1.921fb54442d18p+2,
              -0x1.fc43899911feep+0, -0x1.c7e2cad75cbf8p-1, 0x1.b0ce871368011p+0}, // sheared, whole
        Known{0x1.ec1966a18083bp-4, 0x1.e74d265a2ee3ap+0, -0x1.8933a66b295efp-5,
              0x1.fd9e6a6b4c668p-3, 0x0.0p+0, 0x0.0p+0, -0x1.fe87ec24fa3dap+1, 0x1.6218e7a8a979cp+2,
              -0x1.745e478e85020p+0, -0x1.451f54721e09cp-1, 0x1.68c3596d7fa72p+0}, // sheared, arc
    };

    /**
     *  An ellipse of a row.
     *
     *  @param  known The row.
     *  @return       Its arc.
     */
    auto Arced(const Known &known) -> Ellipse
    {
      return {.Map   = {.X      = {.X = known[0], .Y = known[1]},
                        .Y      = {.X = known[2], .Y = known[3]},
                        .Origin = {.X = known[4], .Y = known[5]}},
              .Start = Radians(known[6]),
              .Sweep = Radians(known[7])};
    }

    /**
     *  A step of a double at the largest of a row's numbers: how near the distance is told.
     *
     *  @param  known The row.
     *  @return       The step.
     */
    auto RowUnit(const Known &known) -> double
    {
      double most = 0.0;
      for (std::size_t index = 0; index < 10; ++index)
      {
        if (index != 6 && index != 7)
        {
          most = std::max(most, std::abs(known.at(index)));
        }
      }
      return 0x1p-52 * most;
    }

    /**
     *  Whether no dense sample of an arc is nearer a point than its nearest place, beyond a few
     *  steps of a double.
     *
     *  @param  ellipse The arc.
     *  @param  point   The point.
     *  @param  unit    A step of a double at the size of the numbers.
     *  @return         Whether none is.
     */
    auto Settled(const Ellipse &ellipse, Point point, double unit) -> bool
    {
      const double  least   = Distance(point, Nearest(ellipse, point).At);
      constexpr int Samples = 20'000;
      for (int index = 0; index <= Samples; ++index)
      {
        if (Distance(point, At(ellipse, static_cast<double>(index) / Samples))
            < least - (4.0 * unit))
        {
          return false;
        }
      }
      return true;
    }
  } // namespace

  TEST_CASE("secular: the nearest place within four steps of a double of mpmath's distance")
  {
    for (std::size_t index = 0; index < Skewed.size(); ++index)
    {
      CAPTURE(index);
      const Known  &known   = Skewed.at(index);
      const Ellipse ellipse = Arced(known);
      const Point   point   = {.X = known[8], .Y = known[9]};
      const Place   place   = Nearest(ellipse, point);
      CHECK(std::abs(Distance(point, place.At) - known[10]) <= 4.0 * RowUnit(known));
      CHECK(place.At == At(ellipse, place.Along));
      CHECK(place.Along >= 0.0);
      CHECK(place.Along <= 1.0);
      CHECK(Settled(ellipse, point, RowUnit(known)));
    }
  }

  TEST_CASE("secular: ties across the major axis, the evolute's cusp, and circles")
  {
    const Place middle = Nearest(Ellipse{.Map = Scaling(3, 1), .Start = {}, .Sweep = Turn}, {});
    CHECK(middle.Along == 0.25);
    CHECK(Distance(middle.At, {.X = 0, .Y = 1}) < 1e-15);

    const Place cusp =
        Nearest(Ellipse{.Map = Scaling(2, 1), .Start = {}, .Sweep = Turn}, {.X = 1.5, .Y = 0});
    CHECK(cusp.Along == 0.0);
    CHECK(cusp.At    == Point{.X = 2, .Y = 0});

    const Ellipse circle = {.Map = Scaling(2), .Start = {}, .Sweep = Turn};
    CHECK(Distance(Nearest(circle, {.X = 3, .Y = 4}).At, {.X = 1.2, .Y = 1.6}) < 1e-15);
    CHECK(Nearest(circle, {}).Along == 0.0);
    CHECK_FALSE(Closest(circle, {.X = 3, .Y = 4}));
    CHECK_FALSE(Closest(Ellipse{.Map = Scaling(3, 1), .Start = {}, .Sweep = {}}, {}));
    CHECK(Nearest(Ellipse{.Map = Scaling(3, 1), .Start = Degrees(90), .Sweep = {}}, {}).At
          == Point{.X = 0, .Y = 1});

    /**
     *  On the evolute's cusp, a² − b² − ax nothing, the start below the root still a cube root.
     */
    const Place onCusp =
        Nearest(Ellipse{.Map = Scaling(2, 1), .Start = {}, .Sweep = Turn}, {.X = 1.5, .Y = 0.01});
    CHECK(onCusp.At.Y > 0.0);
    CHECK(Settled(Ellipse{.Map = Scaling(2, 1), .Start = {}, .Sweep = Turn}, {.X = 1.5, .Y = 0.01},
                  0x1p-52));
  }

  TEST_CASE("secular: an arc that leaves the nearest foot out, by the next least or an end")
  {
    const Place below =
        Nearest(Ellipse{.Map = Scaling(3, 1), .Start = Degrees(180), .Sweep = Degrees(180)},
                {.X = 0.5, .Y = 0.2});
    CHECK(below.At.Y < 0.0);
    CHECK(Distance(below.At, {.X = 0.5, .Y = 0.2})
          == doctest::Approx(1.1838364388632134).epsilon(1e-15));

    const Ellipse most = {.Map = Scaling(3, 1), .Start = {}, .Sweep = Degrees(-300)};
    CHECK(Nearest(most, {.X = 2, .Y = 2}).Along == 1.0);
    const Ellipse almost = {.Map = Scaling(3, 1), .Start = {}, .Sweep = Degrees(-330)};
    CHECK(Nearest(almost, {.X = 2.5, .Y = 0.25}).Along == 1.0);
  }
} // namespace Geoxy
