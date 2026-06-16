from dataclasses import dataclass


@dataclass(frozen=True)
class WindowConfig:
    """Configuration for pre/post-stimulus analysis windows.

    Two modes:
        Train mode (symmetric=False): Standard analysis with pre-stim baseline,
            full stim train, and post-stim window. Used for 700ms, 2700ms.

        Symmetric mode (symmetric=True): Symmetric window around first pulse
            onset. Pre-stim = [-window_ms, 0], stim = [0, +window_ms].
            No post-stim. No blanking gap. Used for control analyses (150ms, etc.)

    Attributes:
        pre_stim_ms: Pre-stimulus baseline window in ms.
        post_stim_ms: Post-stimulus observation window in ms (ignored if symmetric).
        label: Human-readable label used in output filenames.
        clip_to_next_train: If True, clip post-stim window to avoid next train's baseline.
        symmetric: If True, use symmetric window around first pulse onset.
    """
    pre_stim_ms: float = 700
    post_stim_ms: float = 700
    label: str = "700ms"
    clip_to_next_train: bool = True
    symmetric: bool = False

    FS = 30000  # sampling rate in Hz

    @property
    def pre_stim_samples(self) -> int:
        return int(self.pre_stim_ms * self.FS / 1000)

    @property
    def post_stim_samples(self) -> int:
        return int(self.post_stim_ms * self.FS / 1000)

    @classmethod
    def standard_700ms(cls):
        return cls(pre_stim_ms=700, post_stim_ms=700, label="700ms")

    @classmethod
    def extended_2700ms(cls):
        return cls(pre_stim_ms=700, post_stim_ms=2700, label="2700ms")

    @classmethod
    def symmetric_150ms(cls):
        return cls(pre_stim_ms=150, post_stim_ms=0, label="sym150ms",
                   clip_to_next_train=False, symmetric=True)

    @classmethod
    def symmetric_120ms(cls):
        return cls(pre_stim_ms=120, post_stim_ms=0, label="sym120ms",
                   clip_to_next_train=False, symmetric=True)

    @classmethod
    def symmetric_100ms(cls):
        return cls(pre_stim_ms=100, post_stim_ms=0, label="sym100ms",
                   clip_to_next_train=False, symmetric=True)

    @classmethod
    def symmetric_700ms(cls):
        return cls(pre_stim_ms=700, post_stim_ms=0, label="sym700ms",
                   clip_to_next_train=False, symmetric=True)
