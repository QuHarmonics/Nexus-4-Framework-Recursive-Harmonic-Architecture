// Core loop wrapper with observable pressure feedback
using System;
using UniRx;
using UnityEngine;

public class HarmonicLoop : MonoBehaviour
{
    private int maxIterations = 10000;
    private int loopIndex = 0;
    private float harmonicPressure = 9;
    private float resonance = 9;

    // Exposed observable for external systems
    public ReactiveProperty<float> PressureFeedback = new ReactiveProperty<float>();
    public ReactiveProperty<float> ResonanceFeedback = new ReactiveProperty<float>();

    void Start()
    {
        Observable.EveryUpdate()
            .TakeUntilDestroy(this)
            .Subscribe(_ => ExecuteLoop());
    }

    void ExecuteLoop()
    {
        // Simulated internal loop logic (placeholder)
        float fluctuation = Mathf.Sin(loopIndex * 0.05f);
        harmonicPressure = Mathf.Abs(fluctuation);
        resonance = Mathf.PerlinNoise(loopIndex * 0.002f, harmonicPressure);

        // Only update wrapper (not loop logic)
        PressureFeedback.Value = harmonicPressure;
        ResonanceFeedback.Value = resonance;

        loopIndex = (loopIndex + 1) % maxIterations;
    }
}
