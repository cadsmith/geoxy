/***************************************************************************************************
 *  @file      source/radians.hpp
 *  @brief     The constants of angles held in radians: 2 / pi to 1280 bits, pi and its half as
 *             pairs, and the arctangents the headings start from.
 *  @copyright Copyright (c) 2026 EL JABIRI Tarik. Released under the MIT License.
 **************************************************************************************************/

/**
 *  Written with MPFR at 4096 bits.
 */

#ifndef GEOXY_RADIANS_HPP
#define GEOXY_RADIANS_HPP

#include "twofold.hpp"

#include <array>
#include <cstdint>

// NOLINTBEGIN(modernize-use-std-numbers): table entries, not those constants
namespace Geoxy::Table
{
  /**
   *  pi.
   */
  inline constexpr Twofold::Two Pi = {.High = 0x1.921fb54442d18p+1, .Low = 0x1.1a62633145c07p-53};

  /**
   *  pi / 2.
   */
  inline constexpr Twofold::Two HalfPi = {.High = 0x1.921fb54442d18p+0,
                                          .Low  = 0x1.1a62633145c07p-54};

  /**
   *  2 / pi's first 1280 bits after the point, 32 to a word, most significant first.
   */
  inline constexpr std::array<std::uint32_t, 40> TwoOverPi = {
      0xa2f9836eU, 0x4e441529U, 0xfc2757d1U, 0xf534ddc0U, 0xdb629599U, 0x3c439041U, 0xfe5163abU,
      0xdebbc561U, 0xb7246e3aU, 0x424dd2e0U, 0x06492eeaU, 0x09d1921cU, 0xfe1deb1cU, 0xb129a73eU,
      0xe88235f5U, 0x2ebb4484U, 0xe99c7026U, 0xb45f7e41U, 0x3991d639U, 0x835339f4U, 0x9c845f8bU,
      0xbdf9283bU, 0x1ff897ffU, 0xde05980fU, 0xef2f118bU, 0x5a0a6d1fU, 0x6d367ecfU, 0x27cb09b7U,
      0x4f463f66U, 0x9e5fea2dU, 0x7527bac7U, 0xebe5f17bU, 0x3d0739f7U, 0x8a5292eaU, 0x6bfb5fb1U,
      0x1f8d5d08U, 0x56033046U, 0xfc7b6babU, 0xf0cfbc20U, 0x9af4361dU,
  };

  /**
   *  atan(k / 128) in radians, k from 0 to 128.
   */
  inline constexpr std::array<Twofold::Two, 129> ArcRadians = {{
      {.High =               0x0p+0, .Low =                 0x0p+0},
      {.High = 0x1.fffd555bbba97p-8, .Low =  0x1.68062351fbbe6p-63},
      {.High = 0x1.fff555bbb729bp-7, .Low =  -0x1.220c39d4dff5p-61},
      {.High = 0x1.7fee0184a5c36p-6, .Low = -0x1.43189fc0a354bp-60},
      {.High = 0x1.ffd55bba97625p-6, .Low = -0x1.5ec431444912cp-60},
      {.High = 0x1.3fd65f169c9d9p-5, .Low =  0x1.7230a716461b5p-61},
      {.High = 0x1.7fb818430da2ap-5, .Low = -0x1.86ef8f794f105p-63},
      {.High = 0x1.bf8ddf139c444p-5, .Low = -0x1.89fe34b2a7fa8p-59},
      {.High = 0x1.ff55bb72cfdeap-5, .Low = -0x1.c934d86d23f1dp-60},
      {.High = 0x1.1f86dbf082d59p-4, .Low = -0x1.095dc7732ef81p-59},
      {.High = 0x1.3f59f0e7c559dp-4, .Low =  0x1.ac4ce285df847p-58},
      {.High = 0x1.5f2324fd2d7b2p-4, .Low =  0x1.8a8da4401318ep-58},
      {.High = 0x1.7ee182602f10fp-4, .Low = -0x1.cfb654c0c3d98p-58},
      {.High = 0x1.9e94153cfdcf1p-4, .Low =  0x1.a332e1d69c47ep-58},
      {.High = 0x1.be39ebe6f07c3p-4, .Low =  0x1.f7b8f29a05987p-58},
      {.High = 0x1.ddd21701eba6ep-4, .Low =  0x1.94effcd76fe58p-58},
      {.High = 0x1.fd5ba9aac2f6ep-4, .Low = -0x1.cd37686760c17p-59},
      {.High = 0x1.0e6adccf40882p-3, .Low = -0x1.d71a31bb98d0dp-57},
      {.High = 0x1.1e1fafb043727p-3, .Low = -0x1.b485914dacf8cp-59},
      {.High = 0x1.2dcbdb2fba1ffp-3, .Low =  0x1.8f28705561534p-58},
      {.High = 0x1.3d6eee8c6626cp-3, .Low =  0x1.61a3b0ce9281bp-57},
      {.High = 0x1.4d087a9da4f17p-3, .Low =  0x1.1f323f1adf158p-57},
      {.High = 0x1.5c9811e3ec26ap-3, .Low = -0x1.054ab2c010f3dp-58},
      {.High = 0x1.6c1d4898933d9p-3, .Low = -0x1.2954a7603c427p-58},
      {.High = 0x1.7b97b4bce5b02p-3, .Low =  0x1.347b0b4f881cap-58},
      {.High = 0x1.8b06ee2879c29p-3, .Low = -0x1.118cd30308c4fp-57},
      {.High = 0x1.9a6a8e96c8626p-3, .Low =  0x1.cf601e7b4348ep-59},
      {.High = 0x1.a9c231b403279p-3, .Low =  0x1.0e8bbe89cca85p-57},
      {.High = 0x1.b90d7529260a2p-3, .Low =  0x1.17b10d2e0e5abp-61},
      {.High = 0x1.c84bf8a742e6ep-3, .Low = -0x1.95bdd0682ea26p-58},
      {.High = 0x1.d77d5df205736p-3, .Low =  0x1.c648d1534597ep-57},
      {.High = 0x1.e6a148e96ec4dp-3, .Low =  0x1.866b22029f765p-57},
      {.High = 0x1.f5b75f92c80ddp-3, .Low =  0x1.8ab6e3cf7afbdp-57},
      {.High = 0x1.025fa510665b6p-2, .Low = -0x1.672df6832fa48p-56},
      {.High = 0x1.09dc597d86362p-2, .Low =  0x1.62e47390cb865p-56},
      {.High = 0x1.1151a362431cap-2, .Low = -0x1.4dc8dc9077b9fp-56},
      {.High = 0x1.18bf5a30bf178p-2, .Low =  0x1.30ca4748b1bf9p-57},
      {.High = 0x1.2025567e47c96p-2, .Low = -0x1.1832328f4290ep-57},
      {.High = 0x1.278372057ef46p-2, .Low = -0x1.077cdd36dfc81p-56},
      {.High = 0x1.2ed987a823cfep-2, .Low =  0x1.b91258ea012cap-57},
      {.High = 0x1.362773707ebccp-2, .Low = -0x1.963a544b672d8p-57},
      {.High = 0x1.3d6d129271134p-2, .Low =  0x1.137ca41cc958ap-56},
      {.High = 0x1.44aa436c2af0ap-2, .Low = -0x1.5d5e43c55b3bap-56},
      {.High = 0x1.4bdee586890e7p-2, .Low = -0x1.e4dc77c22a757p-57},
      {.High = 0x1.530ad9951cd4ap-2, .Low = -0x1.2566480884082p-57},
      {.High = 0x1.5a2e0175e0f4ep-2, .Low =  0x1.13b7a8f82e457p-56},
      {.High = 0x1.614840309cfe2p-2, .Low =   -0x1.a725715711fp-56},
      {.High = 0x1.685979f5fa6fep-2, .Low = -0x1.257814d1ada9cp-59},
      {.High = 0x1.6f61941e4def1p-2, .Low = -0x1.c63aae6f6e918p-56},
      {.High = 0x1.7660752817502p-2, .Low =   -0x1.dd11791cc76p-59},
      {.High = 0x1.7d5604b63b3f7p-2, .Low =  0x1.69c885c2b249ap-56},
      {.High = 0x1.84422b8df95d7p-2, .Low =  0x1.d76a0299b41b6p-56},
      {.High = 0x1.8b24d394a1b25p-2, .Low =  0x1.b6d0ba3748fa8p-56},
      {.High = 0x1.91fde7cd0c662p-2, .Low =  0x1.1074188054b53p-56},
      {.High = 0x1.98cd5454d6b18p-2, .Low =  0x1.9e6c988fd0a77p-56},
      {.High = 0x1.9f93066168002p-2, .Low = -0x1.c827047c9439ap-56},
      {.High = 0x1.a64eec3cc23fdp-2, .Low = -0x1.24dec1b50b7ffp-56},
      {.High = 0x1.ad00f5422058bp-2, .Low =  0x1.fc4c33891d2e8p-56},
      {.High = 0x1.b3a911da65c6cp-2, .Low =   0x1.ae187b1ca504p-56},
      {.High = 0x1.ba473378624a5p-2, .Low =  0x1.519a1b46e4affp-56},
      {.High =  0x1.c0db4c94ec9fp-2, .Low = -0x1.cc1ce70934c34p-56},
      {.High = 0x1.c76550aad71f9p-2, .Low = -0x1.74b8bff7043e4p-56},
      {.High = 0x1.cde53432c1351p-2, .Low = -0x1.a2cfa4418f1adp-56},
      {.High = 0x1.d45aec9ec862bp-2, .Low =  0x1.89421163ef92dp-57},
      {.High = 0x1.dac670561bb4fp-2, .Low =  0x1.a2b7f222f65e2p-56},
      {.High =  0x1.e127b6b0744bp-2, .Low = -0x1.2b0986398d4abp-58},
      {.High = 0x1.e77eb7f175a34p-2, .Low =  0x1.0e53dc1bf3435p-56},
      {.High = 0x1.edcb6d43f8435p-2, .Low = -0x1.fc976330884e4p-58},
      {.High = 0x1.f40dd0b541418p-2, .Low = -0x1.a3992dc382a23p-57},
      {.High = 0x1.fa45dd3029259p-2, .Low = -0x1.ca563dc28d8b5p-56},
      {.High = 0x1.0039c73c1a40cp-1, .Low = -0x1.b32c949c9d593p-55},
      {.High = 0x1.034b709250488p-1, .Low =   0x1.8f9b38d85541p-56},
      {.High =  0x1.0657e94db30dp-1, .Low = -0x1.d5b495f6349e6p-56},
      {.High =  0x1.095f30861a59p-1, .Low = -0x1.121b20a15a9f3p-56},
      {.High = 0x1.0c6145b5b43dap-1, .Low =  0x1.974fa13b5404fp-58},
      {.High = 0x1.0f5e28b67e295p-1, .Low =   0x1.311b17ec990dp-65},
      {.High = 0x1.1255d9bfbd2a9p-1, .Low = -0x1.2bdaee1c0ee35p-58},
      {.High = 0x1.154859637646ap-1, .Low = -0x1.4ba7c548bf3c3p-55},
      {.High = 0x1.1835a88be7c13p-1, .Low =  0x1.c621cec00c301p-55},
      {.High = 0x1.1b1dc87904285p-1, .Low = -0x1.21e8c8aef8f29p-57},
      {.High = 0x1.1e00babdefeb4p-1, .Low = -0x1.928df287a668fp-58},
      {.High = 0x1.20de813e823b2p-1, .Low = -0x1.791d753ebb744p-55},
      {.High = 0x1.23b71e2cc9e6ap-1, .Low =  0x1.c421c9f38224ep-57},
      {.High = 0x1.268a940696da6p-1, .Low =  0x1.d1348a04c73ccp-58},
      {.High = 0x1.2958e59308e31p-1, .Low = -0x1.09e73b0c6c087p-56},
      {.High = 0x1.2c2215e024466p-1, .Low = -0x1.4b810da3a4be1p-59},
      {.High = 0x1.2ee628406cbcap-1, .Low =  0x1.c5d5e9ff0cf8dp-55},
      {.High = 0x1.31a52048874bep-1, .Low =  0x1.40cab87a7ac24p-55},
      {.High = 0x1.345f01cce37bbp-1, .Low =  0x1.1021137c71102p-55},
      {.High = 0x1.3713d0df6c504p-1, .Low = -0x1.4f789e031606dp-58},
      {.High = 0x1.39c391cd4171ap-1, .Low = -0x1.2304331d8bf46p-55},
      {.High = 0x1.3c6e491c78dc5p-1, .Low = -0x1.e145094fd0ba7p-55},
      {.High = 0x1.3f13fb89e96f4p-1, .Low =   0x1.ecf8b492644fp-56},
      {.High = 0x1.41b4ae06fea41p-1, .Low =  0x1.3d60a53277652p-57},
      {.High = 0x1.445065b795b56p-1, .Low = -0x1.f76d0163f79c8p-56},
      {.High = 0x1.46e727efe4716p-1, .Low = -0x1.39b9b1b844cc9p-57},
      {.High = 0x1.4978fa3269ee1p-1, .Low =  0x1.2419a87f2a458p-56},
      {.High = 0x1.4c05e22de94e5p-1, .Low = -0x1.c0ac1f09f2edfp-55},
      {.High = 0x1.4e8de5bb6ec04p-1, .Low =  0x1.4a33dbeb3796cp-55},
      {.High = 0x1.51110adc5ed81p-1, .Low =  0x1.23dcd6832a63ep-56},
      {.High = 0x1.538f57b89061fp-1, .Low = -0x1.1bb74abda520cp-55},
      {.High = 0x1.5608d29c70c34p-1, .Low =  0x1.9939cf0de8088p-55},
      {.High = 0x1.587d81f732fbbp-1, .Low =  -0x1.5e5c9d8c5a95p-56},
      {.High = 0x1.5aed6c5909517p-1, .Low =  0x1.7312f714a9436p-55},
      {.High = 0x1.5d58987169b18p-1, .Low =  0x1.0028e4bc5e7cap-57},
      {.High = 0x1.5fbf0d0d5cc4ap-1, .Low = -0x1.b4cfd000b7158p-58},
      {.High = 0x1.6220d115d7b8ep-1, .Low = -0x1.2b785350ee8c1p-57},
      {.High =  0x1.647deb8e20b9p-1, .Low = -0x1.eca04023a51cfp-58},
      {.High = 0x1.66d663923e087p-1, .Low = -0x1.6ea6febe8bbbap-56},
      {.High = 0x1.692a40556fb6ap-1, .Low =  0x1.d94b95a8ea2ccp-55},
      {.High = 0x1.6b798920b3d99p-1, .Low = -0x1.a80386188c50ep-55},
      {.High = 0x1.6dc44551553afp-1, .Low = -0x1.bf8863573828ep-58},
      {.High = 0x1.700a7c5784634p-1, .Low = -0x1.8c34d25aadef6p-56},
      {.High = 0x1.724c35b4fae7bp-1, .Low =  0x1.948b32db3499bp-58},
      {.High = 0x1.748978fba8e0fp-1, .Low =  0x1.7b2a6165884a1p-59},
      {.High =  0x1.76c24dcc6c6cp-1, .Low =  0x1.1952551adc83dp-55},
      {.High = 0x1.78f6bbd5d315ep-1, .Low =   0x1.406a08980374p-55},
      {.High = 0x1.7b26cad2e50fep-1, .Low = -0x1.ce80df30411fbp-55},
      {.High = 0x1.7d528289fa093p-1, .Low =  0x1.560821e2f3aa9p-55},
      {.High = 0x1.7f79eacb97898p-1, .Low =  0x1.fd5ca80ead221p-55},
      {.High = 0x1.819d0b7158a4dp-1, .Low = -0x1.bf76229d3b917p-56},
      {.High = 0x1.83bbec5cdee22p-1, .Low =  0x1.3107104ffc6c3p-57},
      {.High = 0x1.85d69576cc2c5p-1, .Low =  0x1.6b66e7fc8b8c3p-57},
      {.High = 0x1.87ed0eadc5a2ap-1, .Low =  0x1.0af5ad957f4bcp-56},
      {.High = 0x1.89ff5ff57f1f8p-1, .Low = -0x1.55b9a5e177a1bp-55},
      {.High = 0x1.8c0d9145cf49dp-1, .Low =  0x1.bea4076dc4333p-55},
      {.High = 0x1.8e17aa99cc05ep-1, .Low = -0x1.ec182ab042f61p-56},
      {.High = 0x1.901db3eeef187p-1, .Low =  0x1.68665e5603c8fp-55},
      {.High = 0x1.921fb54442d18p-1, .Low =  0x1.1a62633145c07p-55},
  }};

  /**
   *  floor(pi 2^4096) in 32-bit words, least significant first.
   */
  inline constexpr std::array<std::uint32_t, 129> PiWords = {
      0x80991b7bU, 0x075372c9U, 0x49f1c09bU, 0xdb0fead3U, 0xd00a1248U, 0x37d0d724U, 0x429b023dU,
      0x1bfedf72U, 0x363f7706U, 0x4ed3aa62U, 0x56c16aa6U, 0x401a449fU, 0x85c12073U, 0xe06f75d8U,
      0x3b8b5ebeU, 0x7d84a5c3U, 0x456f9fb4U, 0x8cee8619U, 0x82430e88U, 0x66ca593eU, 0x39af0176U,
      0xa1f1651dU, 0x7efb2a98U, 0xba3bf050U, 0x137a3be4U, 0x6eef0b6cU, 0xab5133a3U, 0x960fa728U,
      0xd8542f68U, 0x6a51a0d2U, 0xabd388f0U, 0x670c9c61U, 0xf6e96c9aU, 0x21c66842U, 0x9e1f9b5eU,
      0x69c8f04aU, 0xa4842004U, 0x2e0b4482U, 0x83f44239U, 0x0f6d6ff3U, 0xd396acc5U, 0x23893e81U,
      0xeb651b88U, 0xdc262302U, 0xe98575b1U, 0xef845d5dU, 0x5dec8032U, 0x487cac60U, 0xfb21a991U,
      0x61d809ccU, 0x66282193U, 0xc4bfe81bU, 0x6b4bb9afU, 0x3b8f4898U, 0x28958677U, 0x7a325381U,
      0x6c24cf5cU, 0xafd6ba33U, 0x9b87931eU, 0xce5c3e16U, 0x741831f6U, 0x2ba9c55dU, 0x636fbc2aU,
      0xb3ee1411U, 0x7c72e993U, 0xa15486afU, 0x1141e8ceU, 0xb4cc5c34U, 0x2aab10b6U, 0x55ca396aU,
      0x63e81440U, 0x57489862U, 0xaa55ab94U, 0xe65525f3U, 0x55605c60U, 0x78af2fdaU, 0xbd314b27U,
      0xd71577c1U, 0xb01e8a3eU, 0x6c9e0e8bU, 0x603a180eU, 0x8e79dcb0U, 0xb8db38efU, 0xca417918U,
      0x286085f0U, 0xc5d1b023U, 0x2af26013U, 0x9c30d539U, 0xc25a59b5U, 0x7b54a41dU, 0x82154aeeU,
      0x718bcd58U, 0x728eb658U, 0x0d95748fU, 0xf4933d7eU, 0xa458fea3U, 0x71574e69U, 0x636920d8U,
      0x858efc16U, 0x0801f2e2U, 0xb3916cf7U, 0x24a19947U, 0xf12c7f99U, 0xba7c9045U, 0x6a267e96U,
      0xb8e1afedU, 0xd01adfb7U, 0x2ffd72dbU, 0x98dfb5acU, 0xd1310ba6U, 0x8979fb1bU, 0x9216d5d9U,
      0xb5470917U, 0x3f84d5b5U, 0xc97c50ddU, 0xc0ac29b7U, 0x34e90c6cU, 0xbe5466cfU, 0x38d01377U,
      0x452821e6U, 0xec4e6c89U, 0x082efa98U, 0x299f31d0U, 0xa4093822U, 0x03707344U, 0x13198a2eU,
      0x85a308d3U, 0x243f6a88U, 0x00000003U,
  };
} // namespace Geoxy::Table

// NOLINTEND(modernize-use-std-numbers)

#endif /* GEOXY_RADIANS_HPP */
