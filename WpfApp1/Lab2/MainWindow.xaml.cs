/*using System;
using System.Collections.Generic;
using System.Numerics;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Media.Imaging;

namespace Lab2
{
    // --- DATENSTRUKTUREN FÜR RAY TRACING ---

    /// <summary>
    /// Repräsentiert einen Strahl im 3D-Raum.
    /// </summary>
    public struct Ray
    {
        public Vector3 Origin;
        public Vector3 Direction;
    }

    /// <summary>
    /// Definiert die Oberflächeneigenschaften eines Objekts.
    /// Vorerst nur die Farbe (im linearen RGB-Raum).
    /// </summary>
    public struct Material
    {
        public Vector3 Color;
    }

    /// <summary>
    /// Speichert Informationen über einen Ray-Objekt-Schnittpunkt.
    /// </summary>
    public struct HitInfo
    {
        public bool DidHit;
        public float Distance;
        public Vector3 HitPoint;
        public Vector3 Normal;
        public Material Material;
    }

    /// <summary>
    /// Repräsentiert eine Kugel in der 3D-Szene.
    /// </summary>
    public class Sphere
    {
        public Vector3 Center;
        public float Radius;
        public Material Material;

        /// <summary>
        /// Berechnet den Schnittpunkt eines Strahls mit dieser Kugel.
        /// </summary>
        /// <param name="ray">Der Strahl, der getestet wird.</param>
        /// <returns>HitInfo mit den Details des Treffers.</returns>
        public HitInfo Intersect(Ray ray)
        {
            Vector3 oc = ray.Origin - Center;
            float a = Vector3.Dot(ray.Direction, ray.Direction);                         
            float b = 2.0f * Vector3.Dot(oc, ray.Direction);
            float c = Vector3.Dot(oc, oc) - Radius * Radius;
            float discriminant = b * b - 4 * a * c;

            if (discriminant < 0)
            {
                return new HitInfo { DidHit = false }; // Kein Treffer
            }
            
            float t = (-b - MathF.Sqrt(discriminant)) / (2.0f * a);
            // Wir wollen nur Treffer, die vor der Kamera liegen (t > 0)
            if (t < 0.001f) {
                t = (-b + MathF.Sqrt(discriminant)) / (2.0f * a);
                if (t < 0.001f) return new HitInfo { DidHit = false };
            }

            Vector3 hitPoint = ray.Origin + t * ray.Direction;
            return new HitInfo
            {
                DidHit = true,
                Distance = t,
                HitPoint = hitPoint,
                Normal = Vector3.Normalize(hitPoint - Center),
                Material = this.Material
            };
        }
    }

    /// <summary>
    /// Enthält alle Objekte der Szene.
    /// </summary>
    public class Scene
    {
        public List<Sphere> Spheres = new List<Sphere>();
    }


    public partial class MainWindow : Window
    {
        const int W = 512, H = 512;
        private WriteableBitmap _wb;

        public MainWindow()
        {
            InitializeComponent();
            Image canvasImage = new Image();
            this.Content = canvasImage;
            _wb = new WriteableBitmap(W, H, 96, 96, PixelFormats.Bgra32, null);
            canvasImage.Source = _wb;

            // Starte den Rendering-Prozess
            RenderScene();
        }
        
        /// <summary>
        /// Die Hauptfunktion, die die Szene aufbaut und rendert.
        /// </summary>
        void RenderScene()
        {
            // 1. Szene definieren (ähnlich deiner Vorlage)
            // Hinweis: Die Farben sind < 1.0, wie in deinem Bild empfohlen (~80% Intensität)
            var scene = new Scene();
            scene.Spheres.Add(new Sphere { Center = new Vector3(-1001, 0, 0), Radius = 1000, Material = new Material { Color = new Vector3(0.8f, 0.1f, 0.1f) } }); // Linke Wand (Rot)
            scene.Spheres.Add(new Sphere { Center = new Vector3(1001, 0, 0), Radius = 1000, Material = new Material { Color = new Vector3(0.1f, 0.1f, 0.8f) } }); // Rechte Wand (Blau)
            scene.Spheres.Add(new Sphere { Center = new Vector3(0, 1001, 0), Radius = 1000, Material = new Material { Color = new Vector3(0.8f, 0.8f, 0.8f) } }); // Decke (Grau)
            scene.Spheres.Add(new Sphere { Center = new Vector3(0, -1001, 0), Radius = 1000, Material = new Material { Color = new Vector3(0.6f, 0.6f, 0.6f) } }); // Boden (Grau)
            scene.Spheres.Add(new Sphere { Center = new Vector3(0, 0, 1001), Radius = 1000, Material = new Material { Color = new Vector3(0.6f, 0.6f, 0.6f) } }); // Rückwand (Grau)

            scene.Spheres.Add(new Sphere { Center = new Vector3(-0.6f, -0.7f, -0.6f), Radius = 0.3f, Material = new Material { Color = new Vector3(0.9f, 0.9f, 0.1f) } }); // Gelbe Kugel
            scene.Spheres.Add(new Sphere { Center = new Vector3(0.3f, -0.4f, 0.3f), Radius = 0.6f, Material = new Material { Color = new Vector3(0.1f, 0.9f, 0.9f) } }); // Cyan Kugel
            
            // 2. Kamera definieren
            /*Vector3 eye = new Vector3(0, 0, -4);#1#
            Vector3 eye = new Vector3(-0.9f, -0.5f, 0.9f);
            /*Vector3 lookAt = new Vector3(0, 0, 6);#1#
            Vector3 lookAt = new Vector3(0, 0, 0);
            /*float fov = 36.0f;#1#
            float fov = 110.0f;
            Vector3 backgroundColor = new Vector3(0,0,0); // Schwarzer Hintergrund

            int stride = _wb.BackBufferStride;
            byte[] pixels = new byte[stride * H];
            
            // Kamera-Setup für die Strahlenberechnung
            float aspectRatio = (float)W / H;
            float fovRad = MathF.PI / 180.0f * fov;
            float viewportHeight = 2.0f * MathF.Tan(fovRad / 2.0f);
            float viewportWidth = viewportHeight * aspectRatio;

            Vector3 forward = Vector3.Normalize(lookAt - eye);
            Vector3 right = Vector3.Normalize(Vector3.Cross(new Vector3(0, 1, 0), forward));
            Vector3 up = Vector3.Cross(forward, right);

            // 3. Durch jeden Pixel iterieren
            for (int y = 0; y < H; y++)
            {
                for (int x = 0; x < W; x++)
                {
                    // a. Strahl für diesen Pixel erstellen (CreateEyeRay)
                    float u = (x / (float)(W - 1)) * 2.0f - 1.0f; // von -1 bis 1
                    float v = ((H - 1 - y) / (float)(H - 1)) * 2.0f - 1.0f; // von -1 bis 1 (y umkehren)
                    
                    Ray ray = new Ray
                    {
                        Origin = eye,
                        Direction = Vector3.Normalize(forward + u * viewportWidth/2 * right + v * viewportHeight/2 * up)
                    };
                    
                    // b. Szene nach Treffern durchsuchen (FindClosestHitPoint)
                    HitInfo closestHit = new HitInfo { DidHit = false, Distance = float.MaxValue };
                    foreach(var sphere in scene.Spheres)
                    {
                        HitInfo currentHit = sphere.Intersect(ray);
                        if(currentHit.DidHit && currentHit.Distance < closestHit.Distance)
                        {
                            closestHit = currentHit;
                        }
                    }

                    // c. Farbe berechnen (ComputeColor)
                    Vector3 linearColor = closestHit.DidHit ? closestHit.Material.Color : backgroundColor;

                    // d. Pixel in den Puffer schreiben
                    int idx = y * stride + x * 4;
                    WriteBgra(pixels, idx, linearColor);
                }
            }

            _wb.Lock();
            _wb.WritePixels(new Int32Rect(0, 0, W, H), pixels, stride, 0);
            _wb.AddDirtyRect(new Int32Rect(0, 0, W, H));
            _wb.Unlock();
        }

        #region Hilfsfunktionen für Farbe und Pixel
        static byte ToSrgb8(float linear)
        {
            if (linear <= 0f) return 0;
            if (linear >= 1f) return 255;
            float srgb = linear <= 0.0031308f ? 12.92f * linear : 1.055f * MathF.Pow(linear, 1f / 2.4f) - 0.055f;
            return (byte)Math.Clamp((int)MathF.Round(srgb * 255f), 0, 255);
        }

        static void WriteBgra(byte[] buf, int idx, Vector3 linearRgb)
        {
            buf[idx + 0] = ToSrgb8(linearRgb.Z); // B
            buf[idx + 1] = ToSrgb8(linearRgb.Y); // G
            buf[idx + 2] = ToSrgb8(linearRgb.X); // R
            buf[idx + 3] = 255; // A
        }
        #endregion
    }
}*/

using System;
using System.Collections.Generic;
using System.IO;
using System.Numerics;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Media.Imaging;

namespace Lab2
{
    // -------------------- TEXTURE CLASS --------------------
    
    public class Texture
    {
        private Vector3[] _pixels;
        private int _width;
        private int _height;

        public Texture(string filePath)
        {
            LoadFromFile(filePath);
        }

        private void LoadFromFile(string filePath)
        {
            BitmapImage bitmapImage = new BitmapImage(new Uri(filePath, UriKind.RelativeOrAbsolute));
            
            // Convert to a format we can read
            FormatConvertedBitmap convertedBitmap = new FormatConvertedBitmap();
            convertedBitmap.BeginInit();
            convertedBitmap.Source = bitmapImage;
            convertedBitmap.DestinationFormat = PixelFormats.Bgra32;
            convertedBitmap.EndInit();

            _width = convertedBitmap.PixelWidth;
            _height = convertedBitmap.PixelHeight;
            
            int stride = _width * 4;
            byte[] pixelData = new byte[_height * stride];
            convertedBitmap.CopyPixels(pixelData, stride, 0);

            // Convert from sRGB bytes to linear RGB float3
            _pixels = new Vector3[_width * _height];
            for (int y = 0; y < _height; y++)
            {
                for (int x = 0; x < _width; x++)
                {
                    int idx = y * stride + x * 4;
                    byte b = pixelData[idx + 0];
                    byte g = pixelData[idx + 1];
                    byte r = pixelData[idx + 2];
                    
                    // Convert sRGB to linear RGB (gamma correction)
                    _pixels[y * _width + x] = new Vector3(
                        SrgbToLinear(r / 255f),
                        SrgbToLinear(g / 255f),
                        SrgbToLinear(b / 255f)
                    );
                }
            }
        }

        // Convert sRGB to linear RGB
        private static float SrgbToLinear(float srgb)
        {
            if (srgb <= 0.04045f)
                return srgb / 12.92f;
            else
                return MathF.Pow((srgb + 0.055f) / 1.055f, 2.4f);
        }

        // Sample texture with bilinear filtering
        public Vector3 Sample(float u, float v)
        {
            // Wrap UV coordinates
            u = u - MathF.Floor(u);
            v = v - MathF.Floor(v);
            
            // Convert to pixel coordinates
            float x = u * (_width - 1);
            float y = v * (_height - 1);
            
            // Bilinear interpolation
            int x0 = (int)MathF.Floor(x);
            int y0 = (int)MathF.Floor(y);
            int x1 = Math.Min(x0 + 1, _width - 1);
            int y1 = Math.Min(y0 + 1, _height - 1);
            
            float fx = x - x0;
            float fy = y - y0;
            
            Vector3 c00 = _pixels[y0 * _width + x0];
            Vector3 c10 = _pixels[y0 * _width + x1];
            Vector3 c01 = _pixels[y1 * _width + x0];
            Vector3 c11 = _pixels[y1 * _width + x1];
            
            Vector3 c0 = c00 * (1 - fx) + c10 * fx;
            Vector3 c1 = c01 * (1 - fx) + c11 * fx;
            
            return c0 * (1 - fy) + c1 * fy;
        }
    }

    // -------------------- RAY TRACING DATA --------------------

    public struct Ray
    {
        public Vector3 E;
        public Vector3 d;
    }

    public struct Material
    {
        public Vector3 Color;      // albedo/diffuse color (used if Texture is null)
        public Vector3 Emission;   // light emission (for area lights)
        public float Roughness;    // 0=mirror, 1=diffuse
        public Texture Texture;    // optional texture map
    }

    public struct HitInfo
    {
        public bool DidHit;
        public float λ;
        public Vector3 H;
        public Vector3 n;
        public Vector2 UV;         // texture coordinates
        public Material Material;
    }

    public class Sphere
    {
        public Vector3 C;
        public float r;
        public Material Material;

        public HitInfo Intersect(Ray ray)
        {
            Vector3 OC = ray.E - C;
            float a = Vector3.Dot(ray.d, ray.d);
            float b = 2.0f * Vector3.Dot(OC, ray.d);
            float c = Vector3.Dot(OC, OC) - r * r;

            float Δ = b * b - 4.0f * a * c;
            if (Δ < 0.0f)
                return new HitInfo { DidHit = false };

            float sqrtΔ = MathF.Sqrt(Δ);
            float λ0 = (-b - sqrtΔ) / (2.0f * a);
            float λ1 = (-b + sqrtΔ) / (2.0f * a);

            const float EPS = 1e-3f;
            float λhit = λ0;
            if (λhit < EPS)
            {
                λhit = λ1;
                if (λhit < EPS)
                    return new HitInfo { DidHit = false };
            }

            Vector3 H = ray.E + λhit * ray.d;
            Vector3 n = Vector3.Normalize(H - C);

            // Calculate UV coordinates for sphere
            Vector2 uv = CalculateSphereUV(n);

            return new HitInfo
            {
                DidHit = true,
                λ = λhit,
                H = H,
                n = n,
                UV = uv,
                Material = this.Material
            };
        }

        // Calculate UV coordinates from sphere normal
        private Vector2 CalculateSphereUV(Vector3 normal)
        {
            // Spherical coordinates: theta (elevation), phi (azimuth)
            float theta = MathF.Acos(Math.Clamp(normal.Y, -1f, 1f));
            float phi = MathF.Atan2(normal.Z, normal.X);
            
            // Convert to UV (0 to 1 range)
            float u = 0.5f + phi / (2.0f * MathF.PI);
            float v = theta / MathF.PI;
            
            return new Vector2(u, v);
        }
    }

    public class Scene
    {
        public List<Sphere> Spheres = new List<Sphere>();

        public HitInfo Intersect(Ray ray)
        {
            HitInfo closest = new HitInfo { DidHit = false, λ = float.MaxValue };
            foreach (var s in Spheres)
            {
                var hit = s.Intersect(ray);
                if (hit.DidHit && hit.λ < closest.λ)
                    closest = hit;
            }
            return closest;
        }
    }

    // -------------------- PATH TRACING --------------------

    public class PathTracer
    {
        private Random _rng;
        private const int MinBounces = 3;  // Always trace this many
        private const float RRProbability = 0.8f;  // Survival probability

        public PathTracer(int seed = 42)
        {
            _rng = new Random(seed);
        }

        // Main path tracing function with Russian Roulette
        public Vector3 Trace(Ray ray, Scene scene, int depth = 0)
        {
            HitInfo hit = scene.Intersect(ray);
            if (!hit.DidHit)
                return Vector3.Zero; // background

            // Add emission (if surface is a light)
            Vector3 color = hit.Material.Emission;

            // Russian Roulette: probabilistically terminate after MinBounces
            if (depth >= MinBounces)
            {
                float rand = (float)_rng.NextDouble();
                if (rand > RRProbability)
                    return color; // Terminate path
                
                // Continue, but boost by 1/probability to remain unbiased
                // We'll multiply the BRDF result by this compensation factor
            }

            // Diffuse BRDF with importance sampling
            Vector3 diffuse = DiffuseBRDF(ray, hit, scene, depth);
            
            // Apply Russian Roulette compensation
            if (depth >= MinBounces)
                diffuse /= RRProbability;
            
            color += diffuse;

            return color;
        }

        // BRDF with support for both diffuse and specular reflection
        private Vector3 DiffuseBRDF(Ray ray, HitInfo hit, Scene scene, int depth)
        {
            float roughness = hit.Material.Roughness;
            
            // Get the surface color (from texture or solid color)
            Vector3 albedo = hit.Material.Texture != null 
                ? hit.Material.Texture.Sample(hit.UV.X, hit.UV.Y) 
                : hit.Material.Color;
            
            // For very smooth surfaces (low roughness), do specular reflection
            if (roughness < 0.5f)
            {
                // Compute perfect reflection direction: d_r = d - 2(d·n)n
                Vector3 d_r = Vector3.Normalize(ray.d - 2.0f * Vector3.Dot(ray.d, hit.n) * hit.n);
                
                // For glossy materials, perturb the reflection direction based on roughness
                Vector3 sampleDir;
                if (roughness > 0.01f)
                {
                    // Perturb reflection direction for glossy surfaces
                    sampleDir = PerturbDirection(d_r, roughness * roughness);
                }
                else
                {
                    // Perfect mirror reflection
                    sampleDir = d_r;
                }
                
                // Create reflected ray
                Ray scattered = new Ray 
                { 
                    E = hit.H + hit.n * 1e-4f,
                    d = sampleDir 
                };

                // Recursively trace
                Vector3 incomingLight = Trace(scattered, scene, depth + 1);
                
                // For specular: return albedo * incoming light
                // Mix in some diffuse based on roughness
                float specularWeight = 1.0f - roughness * 2.0f; // 1.0 at roughness=0, 0.0 at roughness=0.5
                return albedo * incomingLight * specularWeight;
            }
            else
            {
                // Diffuse material: sample hemisphere with cosine weighting
                Vector3 sampleDir = SampleHemisphereCosine(hit.n);
                
                // Create new ray from hit point
                Ray scattered = new Ray 
                { 
                    E = hit.H + hit.n * 1e-4f,
                    d = sampleDir 
                };

                // Recursively trace
                Vector3 incomingLight = Trace(scattered, scene, depth + 1);

                // BRDF for Lambertian: ρ/π where ρ is albedo
                // Monte Carlo estimator with cosine-weighted sampling:
                // The pdf cancels with the cosine term, leaving just ρ * L
                
                return albedo * incomingLight;
            }
        }
        
        // Perturb a direction for glossy reflections
        private Vector3 PerturbDirection(Vector3 direction, float roughness)
        {
            // Sample a point on the hemisphere around the direction
            float r1 = (float)_rng.NextDouble();
            float r2 = (float)_rng.NextDouble();

            // Use roughness to control the cone angle
            float phi = 2.0f * MathF.PI * r1;
            float cosTheta = MathF.Pow(1.0f - r2, 1.0f / (1.0f + roughness * 10.0f));
            float sinTheta = MathF.Sqrt(1.0f - cosTheta * cosTheta);

            // Local coordinates around the reflection direction
            Vector3 tangent = GetPerpendicular(direction);
            Vector3 bitangent = Vector3.Cross(direction, tangent);

            Vector3 perturbed = 
                tangent * (MathF.Cos(phi) * sinTheta) +
                bitangent * (MathF.Sin(phi) * sinTheta) +
                direction * cosTheta;

            return Vector3.Normalize(perturbed);
        }

        // Sample direction in hemisphere with cosine weighting
        private Vector3 SampleHemisphereCosine(Vector3 normal)
        {
            // Uniformly sample unit disk, then project to hemisphere
            float r1 = (float)_rng.NextDouble();
            float r2 = (float)_rng.NextDouble();

            float phi = 2.0f * MathF.PI * r1;
            float sinTheta = MathF.Sqrt(r2);
            float cosTheta = MathF.Sqrt(1.0f - r2);

            // Local coordinates
            Vector3 tangent = GetPerpendicular(normal);
            Vector3 bitangent = Vector3.Cross(normal, tangent);

            Vector3 direction = 
                tangent * (MathF.Cos(phi) * sinTheta) +
                bitangent * (MathF.Sin(phi) * sinTheta) +
                normal * cosTheta;

            return Vector3.Normalize(direction);
        }

        // Get any perpendicular vector to n
        private Vector3 GetPerpendicular(Vector3 n)
        {
            Vector3 a = MathF.Abs(n.X) > 0.9f 
                ? new Vector3(0, 1, 0) 
                : new Vector3(1, 0, 0);
            return Vector3.Normalize(Vector3.Cross(n, a));
        }
    }

    // -------------------- WPF WINDOW / RENDERER --------------------

    public partial class MainWindow : Window
    {
        const int W = 512, H = 512;
        const int SamplesPerPixel = 100; // Increase for better quality
        private WriteableBitmap _wb;

        public MainWindow()
        {
            InitializeComponent();
            var img = new Image();
            Content = img;
            _wb = new WriteableBitmap(W, H, 96, 96, PixelFormats.Bgra32, null);
            img.Source = _wb;

            RenderScene();
        }

        private static (Vector3 f, Vector3 r, Vector3 u, float vw, float vh) 
            BuildCameraBasis(Vector3 E, Vector3 lookAt, float fovDegrees, float aspect)
        {
            Vector3 f = Vector3.Normalize(lookAt - E);
            Vector3 upWorld = new Vector3(0, 1, 0);
            if (MathF.Abs(Vector3.Dot(f, upWorld)) > 0.999f) 
                upWorld = new Vector3(0, 0, 1);
            
            Vector3 r = Vector3.Normalize(Vector3.Cross(upWorld, f));
            Vector3 u = Vector3.Cross(f, r);

            float α = MathF.PI * fovDegrees / 180.0f;
            float vh = 2.0f * MathF.Tan(α * 0.5f);
            float vw = vh * aspect;

            return (f, r, u, vw, vh);
        }

        private static Ray CreateEyeRay(int x, int y, int W, int H,
            Vector3 E, (Vector3 f, Vector3 r, Vector3 u, float vw, float vh) cam)
        {
            float xNorm = (x / (float)(W - 1)) * 2.0f - 1.0f;
            float yNorm = ((H - 1 - y) / (float)(H - 1)) * 2.0f - 1.0f;

            float β = xNorm * (cam.vw * 0.5f);
            float ω = yNorm * (cam.vh * 0.5f);

            Vector3 d = cam.f + β * cam.r + ω * cam.u;
            d = Vector3.Normalize(d);

            return new Ray { E = E, d = d };
        }

        void RenderScene()
        {
            var scene = new Scene();
    
            // ============== YOUR ORIGINAL CORNELL BOX SCENE ==============
            /*
            // Cornell box walls (diffuse)
            scene.Spheres.Add(new Sphere { 
                C = new Vector3(-1001, 0, 0), r = 1000, 
                Material = new Material { 
                    Color = new Vector3(0.8f, 0.1f, 0.1f), 
                    Emission = Vector3.Zero,
                    Roughness = 1.0f,  // fully diffuse
                    Texture = null
                }
            });
            scene.Spheres.Add(new Sphere { 
                C = new Vector3(1001, 0, 0), r = 1000, 
                Material = new Material { 
                    Color = new Vector3(0.1f, 0.1f, 0.8f), 
                    Emission = Vector3.Zero,
                    Roughness = 1.0f,  // fully diffuse
                    Texture = null
                }
            });
    
            // CEILING AS LIGHT SOURCE
            scene.Spheres.Add(new Sphere { 
                C = new Vector3(0, 1001, 0), r = 1000, 
                Material = new Material { 
                    Color = new Vector3(0.8f, 0.8f, 0.8f), 
                    Emission = new Vector3(1.5f, 1.5f, 1.5f), // emissive ceiling
                    Roughness = 1.0f,  // fully diffuse
                    Texture = null
                }
            });
    
            scene.Spheres.Add(new Sphere { 
                C = new Vector3(0, -1001, 0), r = 1000, 
                Material = new Material { 
                    Color = new Vector3(0.6f, 0.6f, 0.6f), 
                    Emission = Vector3.Zero,
                    Roughness = 1.0f,  // fully diffuse
                    Texture = null
                }
            });
            scene.Spheres.Add(new Sphere { 
                C = new Vector3(0, 0, 1001), r = 1000, 
                Material = new Material { 
                    Color = new Vector3(0.6f, 0.6f, 0.6f), 
                    Emission = Vector3.Zero,
                    Roughness = 1.0f,  // fully diffuse
                    Texture = null
                }
            });

            // Objects - yellow sphere is mirror-like, cyan is glossy
            scene.Spheres.Add(new Sphere { 
                C = new Vector3(-0.6f, -0.7f, -0.6f), r = 0.3f, 
                Material = new Material { 
                    Color = new Vector3(0.9f, 0.9f, 0.1f), 
                    Emission = Vector3.Zero,
                    Roughness = 0.0f,  // mirror (perfect reflection)
                    Texture = null
                }
            });
            scene.Spheres.Add(new Sphere { 
                C = new Vector3(0.3f, -0.4f, 0.3f), r = 0.6f, 
                Material = new Material { 
                    Color = new Vector3(1f, 1f, 1f),
                    Emission = Vector3.Zero,
                    Roughness = 1.0f,  // glossy/partially reflective
                    Texture = null
                }
            });

            // Camera (updated to match your settings)
            Vector3 E = new Vector3(-0.0f, -0.0f, -2f);
            Vector3 lookAt = new Vector3(0, 0, 0);
            float fov = 90f;
            Vector3 backgroundColor = new Vector3(0, 0, 0);
            */
            
            // ============== NEW TEXTURED SCENE DEMO ==============
            
            // Simple white backdrop/floor
            scene.Spheres.Add(new Sphere { 
                C = new Vector3(0, -100.5f, 0), r = 100f, 
                Material = new Material { 
                    Color = new Vector3(0.8f, 0.8f, 0.8f), 
                    Emission = Vector3.Zero,
                    Roughness = 1.0f,  // fully diffuse
                    Texture = null
                }
            });

            // Create a grid of spheres with different materials and roughness values
            // Row 1: Varying roughness (left to right: 0.0 to 1.0)
            for (int i = 0; i < 5; i++)
            {
                float x = -2f + i * 1.0f;
                float roughness = i / 4.0f; // 0.0, 0.25, 0.5, 0.75, 1.0
                
                // Hue varies across the row
                float hue = i / 4.0f;
                Vector3 color = HsvToRgb(hue * 360f, 0.8f, 0.8f);
                
                scene.Spheres.Add(new Sphere { 
                    C = new Vector3(x, 0.5f, 2f), 
                    r = 0.4f, 
                    Material = new Material { 
                        Color = color,
                        Emission = Vector3.Zero,
                        Roughness = roughness,
                        Texture = null
                    }
                });
            }

            // ===== TEXTURED SPHERES - Uncomment and adjust paths when you have textures =====
            
            // Use paths relative to the executable (which is in bin/Debug/net9.0-windows/)
            // Go back up to the Lab2 directory
            string basePath = System.IO.Path.Combine(AppDomain.CurrentDomain.BaseDirectory, @"..\..\..\");
            Texture earthTexture = new Texture(System.IO.Path.Combine(basePath, "coast_sand_rocks_02_diff_1k.jpg"));
            Texture woodTexture = new Texture(System.IO.Path.Combine(basePath, "dark_wood_diff_1k.jpg"));
            Texture checkerTexture = new Texture(System.IO.Path.Combine(basePath, "floor_tiles_06_diff_1k.jpg"));
            
            // Example 1: Diffuse textured sphere
            scene.Spheres.Add(new Sphere { 
                C = new Vector3(-1f, 0.5f, 0f), 
                r = 0.4f, 
                Material = new Material { 
                    Color = new Vector3(1, 1, 1), // Will be multiplied with texture
                    Emission = Vector3.Zero,
                    Roughness = 1.0f,  // diffuse with texture
                    Texture = earthTexture
                }
            });
            
            // Example 2: Glossy textured sphere
            scene.Spheres.Add(new Sphere { 
                C = new Vector3(0f, 0.5f, 0f), 
                r = 0.4f, 
                Material = new Material { 
                    Color = new Vector3(1, 1, 1),
                    Emission = Vector3.Zero,
                    Roughness = 0.9f,  // glossy with texture
                    Texture = woodTexture
                }
            });
            
            // Example 3: Mirror-like textured sphere
            scene.Spheres.Add(new Sphere { 
                C = new Vector3(1f, 0.5f, 0f), 
                r = 0.4f, 
                Material = new Material { 
                    Color = new Vector3(1, 1, 1),
                    Emission = Vector3.Zero,
                    Roughness = 0.9f,  // mirror-like with texture
                    Texture = checkerTexture
                }
            });
           

            // Light dome - large sphere surrounding the entire scene
            scene.Spheres.Add(new Sphere { 
                C = new Vector3(0, 0, 0),  // Center at origin
                r = 200f,  // Large radius to encompass everything
                Material = new Material { 
                    Color = new Vector3(1f, 1f, 1f), 
                    Emission = new Vector3(.5f, .5f, .5f),  // Bright dome light
                    Roughness = 1.0f,
                    Texture = null
                }
            });

            // Camera setup
            Vector3 E = new Vector3(0f, 1.5f, -4f);  // Position camera to view the grid
            Vector3 lookAt = new Vector3(0, 0.5f, 2f);
            float fov = 60f;
            Vector3 backgroundColor = new Vector3(0.1f, 0.1f, 0.15f);  // Dark blue background
            
            float aspect = (float)W / H;
            var cam = BuildCameraBasis(E, lookAt, fov, aspect);

            var tracer = new PathTracer();
            int stride = _wb.BackBufferStride;
            byte[] pixels = new byte[stride * H];

            // Render with multiple samples per pixel
            for (int y = 0; y < H; y++)
            {
                for (int x = 0; x < W; x++)
                {
                    Vector3 colorAccum = Vector3.Zero;
                    
                    for (int s = 0; s < SamplesPerPixel; s++)
                    {
                        Ray ray = CreateEyeRay(x, y, W, H, E, cam);
                        colorAccum += tracer.Trace(ray, scene);
                    }

                    // Average samples
                    Vector3 color = colorAccum / SamplesPerPixel;

                    int idx = y * stride + x * 4;
                    WriteBGRA(pixels, idx, color);
                }

                // Progress update
                if (y % 50 == 0)
                {
                    Title = $"Rendering... {100 * y / H}%";
                }
            }

            _wb.Lock();
            _wb.WritePixels(new Int32Rect(0, 0, W, H), pixels, stride, 0);
            _wb.AddDirtyRect(new Int32Rect(0, 0, W, H));
            _wb.Unlock();

            Title = "Render Complete!";
        }

        static byte ToSrgb8(float linear)
        {
            if (linear <= 0f) return 0;
            if (linear >= 1f) return 255;
            float srgb = linear <= 0.0031308f
                ? 12.92f * linear
                : 1.055f * MathF.Pow(linear, 1f / 2.4f) - 0.055f;
            return (byte)Math.Clamp((int)MathF.Round(srgb * 255f), 0, 255);
        }

        static void WriteBGRA(byte[] buf, int idx, Vector3 linearRgb)
        {
            buf[idx + 0] = ToSrgb8(linearRgb.Z);
            buf[idx + 1] = ToSrgb8(linearRgb.Y);
            buf[idx + 2] = ToSrgb8(linearRgb.X);
            buf[idx + 3] = 255;
        }

        // Convert HSV to RGB (returns linear RGB)
        static Vector3 HsvToRgb(float h, float s, float v)
        {
            float c = v * s;
            float x = c * (1 - MathF.Abs((h / 60f) % 2 - 1));
            float m = v - c;

            float r = 0, g = 0, b = 0;
            if (h < 60) { r = c; g = x; b = 0; }
            else if (h < 120) { r = x; g = c; b = 0; }
            else if (h < 180) { r = 0; g = c; b = x; }
            else if (h < 240) { r = 0; g = x; b = c; }
            else if (h < 300) { r = x; g = 0; b = c; }
            else { r = c; g = 0; b = x; }

            return new Vector3(r + m, g + m, b + m);
        }
    }
}
